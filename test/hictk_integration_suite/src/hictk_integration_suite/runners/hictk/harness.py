# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from datetime import timedelta
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from hictk_integration_suite.runners import Runner


class HictkTestHarness:
    def __init__(
        self,
        hictk_exec: pathlib.Path | str,
        cwd: pathlib.Path | str | None = None,
        tmpdir: pathlib.Path | str | None = None,
    ):
        self._exec = pathlib.Path(hictk_exec) if hictk_exec else None
        self._cwd = pathlib.Path(cwd) if cwd else None
        self._tmpdir = pathlib.Path(tmpdir) if tmpdir else None

        self._id = None
        self._args = []
        self._expect_failure = None

        self._returncode = None
        self._stdout = None
        self._stderr = None

        self._failures = {}
        self._title = None
        self._duration = None

    def _get_hictk_keyword_option(self, option_name: str, default=None) -> Any:
        try:
            i = self._args.index(option_name)
        except ValueError:
            return default

        if i + 1 == len(self._args):
            return default

        return self._args[i + 1]

    def _get_hictk_flag_value(self, flag_name: str) -> bool:
        return flag_name in self._args

    def _run_hictk(
        self,
        args_: List[str],
        timeout: int = 1,
        env_variables: Dict[str, str] | None = None,
        colnames: List[str] | str | None = None,
    ):
        with Runner(self._exec, args_, cwd=self._cwd, tmpdir=self._tmpdir) as runner:
            self._args = runner.args
            self._returncode, self._stdout, self._stderr = runner.run(
                timeout=timeout, env_variables=env_variables, colnames=colnames
            )

    def _validate(self, expect_failure: bool):
        raise NotImplementedError

    def run(
        self,
        args: List[str],
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        expect_failure: bool = False,
        title: str | None = None,
        id: str | None = None,  # noqa
    ) -> Dict[str, Any]:
        if title is None:
            title = str(self)

        self.clear()
        self._id = id
        self._title = title
        self._args = args
        self._expect_failure = expect_failure

        t0 = timer()
        self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        self._validate(expect_failure=expect_failure)
        self._duration = timer() - t0

        return self.status()

    def ok(self) -> bool:
        return self._returncode == 0 and len(self._failures) == 0

    @property
    def returncode(self) -> int:
        return self._returncode

    @property
    def duration(self) -> float:
        return self._duration

    def stderr(self, max_length: int | None = None) -> str:
        payload = "".join(self._stderr)
        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    def stdout(self, max_length: int | None = None) -> str:
        if isinstance(self._stdout, pd.DataFrame):
            if len(self._stdout) == 0:
                payload = ""
            else:
                payload = str(self._stdout)
        else:
            payload = "".join(self._stdout)

        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    @property
    def args(self) -> List[str]:
        return self._args

    @property
    def failures(self) -> Dict[str, str]:
        return self._failures

    def clear(self):
        self._id = None
        self._args = []
        self._expect_failure = None
        self._returncode = None
        self._stdout = None
        self._stderr = None
        self._failures = {}
        self._title = None

    def status(self) -> Dict[str, Any]:
        s = {
            "id": str(self._id),
            "title": str(self._title),
            "args": self.args[1:],
            "elapsed-time": str(timedelta(seconds=self._duration)),
            "exit-code": self._returncode,
            "expect-failure": self._expect_failure,
            "notes": [],
            "status": "PASS" if len(self._failures) == 0 else "FAIL",
        }

        if s["status"] == "PASS":
            return s

        for k, v in self._failures.items():
            if len(v) == 0:
                s["notes"].append(k)
            else:
                s["notes"].append(f"{k}: {v}")

        return s
