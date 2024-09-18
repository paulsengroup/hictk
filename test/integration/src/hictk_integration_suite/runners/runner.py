# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
import shutil
import subprocess as sp
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Tuple

import pandas as pd


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
        logging.debug(f'removing temporary folder "{self._tmpdir}"...')
        os.rmdir(self._tmpdir)

    def _cmd_args(self, include_exec: bool = True) -> List[str]:
        if include_exec:
            args_ = [str(self._exec)]
        else:
            args_ = []
        return args_ + [str(v).strip() for v in self._args]

    @staticmethod
    def _read_table(handle, names: List[str] | None = None) -> pd.DataFrame | str:
        if handle is None:
            raise ValueError("stream cannot be None")

        if len(names) == 0:
            raise ValueError("names cannot be an empty list")

        try:
            logging.debug("reading table from FIFO...")
            df = pd.read_table(handle, names=names)
            logging.debug(f"read {len(df)} records from FIFO")
        except (pd.errors.ParserError, ValueError, OSError) as e:
            logging.warning(f"failed to read table from FIFO: {e}")
            return str(e)

        return df

    @staticmethod
    def _read_file(handle) -> List[str]:
        data = []
        if handle is None:
            return data

        logging.debug("reading data from FIFO...")
        try:
            for line in handle:
                data.append(line)
            logging.debug(f"read {len(data)} records from FIFO...")
        except (ValueError, OSError) as e:
            logging.warning(f"failed to read data from FIFO: {e}")
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
        logging.debug(f"launching subprocess {self._cmd_args()}...")
        return sp.Popen(
            self._cmd_args(),
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            encoding=encoding,
            env=env_variables,
            cwd=self._cwd,
        )

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
    ) -> Tuple[int, str, str]:

        if timeout <= 0:
            raise ValueError("timeout should be a positive number")

        if env_variables is None:
            env_variables = os.environ.copy()

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

            logging.debug(f"waiting for subprocess to return for up to {timeout}s...")
            returncode = proc.wait(timeout)
            logging.debug(f"subprocess terminated with exit code {returncode}")

            delta = time.time() - t0
            time_left = max(0.0, timeout - delta)
            logging.debug(f"waiting for stdout and stderr parsers to return for up to {time_left:.0f}s...")
            return returncode, stdout.result(time_left), stderr.result(time_left)