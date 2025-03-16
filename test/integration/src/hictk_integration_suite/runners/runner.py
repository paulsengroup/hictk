# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
import shutil
import subprocess as sp
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Set, Tuple

import pandas as pd
import structlog


class Runner:
    def __init__(
        self,
        exec: pathlib.Path,
        args_: List,
        cwd: str | None = None,
        tmpdir: pathlib.Path | None = None,
    ):
        self._tmpdir = tempfile.mkdtemp(dir=tmpdir)
        self._exec = shutil.which(exec)
        self._cwd = cwd
        self._args = [str(x) for x in args_]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        structlog.get_logger().debug(f'removing temporary folder "{self._tmpdir}"...')
        os.rmdir(self._tmpdir)

    def _cmd_args(self, include_exec: bool = True) -> List[str]:
        if include_exec:
            args_ = [str(self._exec)]
        else:
            args_ = []
        return args_ + [str(v).strip() for v in self._args]

    @staticmethod
    def _collect_paths(path: pathlib.Path | str | None) -> Set[pathlib.Path]:
        if path is None:
            return set()
        path = pathlib.Path(path)
        paths = set()
        for file in path.iterdir():
            if file.is_dir():
                paths |= Runner._collect_paths(file)
            else:
                paths.add(file)

        return paths

    @staticmethod
    def _generate_path_whitelist(path: pathlib.Path | str | None) -> Set[pathlib.Path]:
        if path is None:
            return set()

        path = pathlib.Path(path)
        structlog.get_logger().debug(f'collecting files under folder "{path}"')

        assert path.is_dir()
        return Runner._collect_paths(path)

    @staticmethod
    def _clean_folder(path: pathlib.Path | str | None, whitelist: Set[pathlib.Path] | None):
        if path is None:
            return

        if whitelist is None:
            whitelist = set()

        path = pathlib.Path(path)
        structlog.get_logger().debug(f'cleaning folder "{path}"')

        assert path.is_dir()
        paths = Runner._collect_paths(path)
        new_paths = paths ^ whitelist
        for path in new_paths:
            path.unlink()

    @staticmethod
    def _read_table(handle, names: List[str] | None = None) -> pd.DataFrame | str:
        if handle is None:
            raise ValueError("stream cannot be None")

        if names is not None and len(names) == 0:
            raise ValueError("names cannot be an empty list")

        try:
            structlog.get_logger().debug("reading table from FIFO...")
            df = pd.read_table(handle, names=names)
            structlog.get_logger().debug(f"read {len(df)} records from FIFO")
        except (pd.errors.ParserError, ValueError, OSError) as e:
            structlog.get_logger().warning(f"failed to read table from FIFO: {e}")
            return str(e)

        return df

    @staticmethod
    def _read_file(handle) -> List[str]:
        data = []
        if handle is None:
            return data

        structlog.get_logger().debug("reading data from FIFO...")
        try:
            for line in handle:
                data.append(line)
            structlog.get_logger().debug(f"read {len(data)} records from FIFO...")
        except (ValueError, OSError) as e:
            structlog.get_logger().warning(f"failed to read data from FIFO: {e}")
            return data

        return data

    def _launch_subproc(
        self,
        stdin,
        stdout,
        stderr,
        encoding: str,
        env_variables: Dict[str, str],
    ) -> sp.Popen:
        structlog.get_logger().debug(f"launching subprocess {self._cmd_args()}...")
        return sp.Popen(
            self._cmd_args(),
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            encoding=encoding,
            env=env_variables,
            cwd=self._cwd,
        )

    def _log_attempt_failure(self, stdout, stderr, timeout: float, attempt_num: int, max_attempts: int):
        structlog.get_logger().warning(
            f"subprocess failed to complete within {timeout} seconds (attempt {attempt_num}/{max_attempts})"
        )

        def log_output(data, label):
            if isinstance(data, pd.DataFrame):
                structlog.get_logger().warning(f"read {len(data)} lines from {label}")
                return

            if len(data) == 0:
                structlog.get_logger().warning(f"{label}: read no data")
                return

            if not isinstance(data, list):
                data = [data]

            i = min(500, len(data))
            for line in data[-i:]:
                structlog.get_logger().warning(f"{label}: {line.strip()}")

        log_output(stdout.result(5.0), "stdout")
        log_output(stderr.result(5.0), "stderr")

    @property
    def args(self):
        return self._cmd_args()

    def run(
        self,
        timeout: float,
        stdin=sp.DEVNULL,
        encoding: str = "utf-8",
        colnames: List[str] | str | None = None,
        env_variables: Dict[str, str] | None = None,
        max_attempts: int = 1,
    ) -> Tuple[int, str, str]:

        if timeout <= 0:
            raise ValueError("timeout should be a positive number")

        if env_variables is None:
            env_variables = os.environ.copy()

        path_whitelist = self._generate_path_whitelist(self._cwd)
        for attempt in range(1, max_attempts + 1):
            self._clean_folder(self._cwd, path_whitelist)
            t0 = time.time()
            with (
                ThreadPoolExecutor(2) as tpool,
                self._launch_subproc(stdin, sp.PIPE, sp.PIPE, encoding, env_variables) as proc,
            ):
                if colnames is None:
                    stdout = tpool.submit(self._read_file, proc.stdout)
                else:
                    if isinstance(colnames, str):
                        assert colnames == "infer"
                        colnames = None

                    stdout = tpool.submit(self._read_table, proc.stdout, colnames)

                stderr = tpool.submit(self._read_file, proc.stderr)

                try:
                    structlog.get_logger().debug(f"waiting for subprocess to return for up to {timeout}s...")
                    returncode = proc.wait(timeout)
                    structlog.get_logger().debug(f"subprocess terminated with exit code {returncode}")
                except sp.TimeoutExpired:
                    proc.kill()
                    self._log_attempt_failure(stdout, stderr, timeout, attempt, max_attempts)
                    if attempt < max_attempts:
                        structlog.get_logger().warning(
                            f"attempting to run subprocess one more time ({max_attempts - attempt} attempt(s) left)"
                        )
                        continue

                    raise

                delta = time.time() - t0
                time_left = max(30.0, timeout - delta)
                structlog.get_logger().debug(
                    f"waiting for stdout and stderr parsers to return for up to {time_left:.0f}s..."
                )
                stdout = stdout.result(time_left)
                stderr = stderr.result(1)
                return returncode, stdout, stderr
