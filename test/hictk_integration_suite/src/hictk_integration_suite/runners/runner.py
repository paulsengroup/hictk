# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import contextlib
import logging
import multiprocessing as mp
import os
import pathlib
import shutil
import stat
import subprocess as sp
import tempfile
import time
from typing import Dict, List, Tuple

import pandas as pd


class Fifo:
    def __init__(self, path: pathlib.Path | str):
        if not stat.S_ISFIFO(os.stat(path).st_mode):
            raise ValueError("path does not point to an existing FIFO")

        self.name = pathlib.Path(path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.remove(self.name)

    def write_handle(self):
        return open(self.name, "a")

    def read_handle(self):
        return open(self.name, "r")

    def name(self) -> pathlib.Path:
        return self.name


class Runner:
    def __init__(self, exec_: pathlib.Path, args: List, cwd: str | None = None, tmpdir: pathlib.Path | None = None):
        self.tmpdir = tempfile.mkdtemp(dir=tmpdir)
        self.exec = shutil.which(exec_)
        self.cwd = cwd
        self.args_ = args

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        logging.debug(f'removing temporary folder "{self.tmpdir}"...')
        os.rmdir(self.tmpdir)

    def _mkfifo(self) -> Fifo:
        with tempfile.NamedTemporaryFile(dir=self.tmpdir, delete_on_close=True) as tmpfile:
            path = tmpfile.name

        logging.debug(f'creating FIFO "{path}"...')
        os.mkfifo(path)
        return Fifo(path)

    def _cmd_args(self, include_exec: bool = True) -> List[str]:
        if include_exec:
            args = [str(self.exec)]
        else:
            args = []
        return args + [str(v).strip() for v in self.args_]

    @staticmethod
    def _read_table(f, names: List[str] | None = None) -> pd.DataFrame | str:
        if f is None:
            raise ValueError("stream cannot be None")

        if len(names) == 0:
            raise ValueError("names cannot be an empty list")

        logging.debug(f'reading table from file "{f}"...')
        try:
            df = pd.read_table(f, names=names)
            logging.debug(f'read {len(df)} records from "{f}"')
        except pd.errors.ParserError as e:
            logging.warning(f'failed to read table from file "{f}: {e}')
            return str(e)

        return df

    @staticmethod
    def _read_file(f: pathlib.Path | None) -> List[str]:
        logging.debug(f'reading data from file "{f}"...')
        if f is None:
            data = []
        else:
            with open(f, "r") as h:
                data = h.readlines()
        logging.debug(f'read {len(data)} records from file "{f}"...')

        return data

    def _launch_subproc(
        self,
        stdin,
        stdout,
        stderr,
        encoding: str,
        env_variables: Dict[str, str],
        ctx: contextlib.ExitStack | None = None,
    ) -> sp.Popen:
        logging.debug(f"launching subprocess {self._cmd_args()}...")
        proc = sp.Popen(
            self._cmd_args(),
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            encoding=encoding,
            env=env_variables,
        )
        if ctx:
            ctx.enter_context(proc)
        return proc

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
        with contextlib.ExitStack() as stack:
            ppool = mp.Pool(2)
            stack.enter_context(ppool)
            stdout_fifo = self._mkfifo()
            stack.enter_context(stdout_fifo)
            stderr_fifo = self._mkfifo()
            stack.enter_context(stderr_fifo)

            if colnames is None:
                stdout_async = ppool.apply_async(self._read_file, args=(stdout_fifo.name,))
            else:
                if isinstance(colnames, str):
                    assert colnames == "infer"
                    colnames = None
                stdout_async = ppool.apply_async(self._read_table, args=(stdout_fifo.name, colnames))
            stderr_async = ppool.apply_async(self._read_file, args=(stderr_fifo.name,))

            with stdout_fifo.write_handle() as stdout, stderr_fifo.write_handle() as stderr:
                proc = self._launch_subproc(stdin, stdout, stderr, encoding, env_variables, ctx=stack)
                logging.debug(f"waiting for subprocess to return for up to {timeout}s...")
                returncode = proc.wait(timeout)
                logging.debug(f"subprocess terminated with exit code {returncode}")

            logging.debug(f"process returned exit code {returncode}")
            delta = time.time() - t0
            time_left = max(0.0, timeout - delta)
            logging.debug(f"waiting for stdout and stderr parsers to return for up to {time_left:.0f}s...")
            return returncode, stdout_async.get(time_left), stderr_async.get(time_left)
