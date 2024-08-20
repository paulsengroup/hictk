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
        self.exec = pathlib.Path(hictk_exec) if hictk_exec else None
        self.cwd = pathlib.Path(cwd) if cwd else None
        self.tmpdir = pathlib.Path(tmpdir) if tmpdir else None

        self.args_ = []

        self.returncode_ = None
        self.stdout_ = None
        self.stderr_ = None

        self.failures_ = {}
        self.title_ = None
        self.duration_ = None

    def _get_hictk_keyword_option(self, option_name: str, default=None) -> Any:
        try:
            i = self.args_.index(option_name)
        except ValueError:
            return default

        if i + 1 == len(self.args_):
            return default

        return self.args_[i + 1]

    def _get_hictk_flag_value(self, flag_name: str) -> bool:
        return flag_name in self.args_

    def _run_hictk(
        self,
        args: List[str],
        timeout: int = 1,
        env_variables: Dict[str, str] | None = None,
        colnames: List[str] | str | None = None,
    ):
        runner = Runner(self.exec, args)
        self.args_ = runner.args()
        self.returncode_, self.stdout_, self.stderr_ = runner.run(
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
    ) -> Dict[str, Any]:
        if title is None:
            title = str(self)

        self.clear()
        self.title_ = title
        self.args_ = args

        t0 = timer()
        self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        self._validate(expect_failure=expect_failure)
        self.duration_ = timer() - t0

        return self.status()

    def ok(self) -> bool:
        return self.returncode_ == 0 and len(self.failures_) == 0

    def returncode(self) -> int:
        return self.returncode_

    def duration(self) -> float:
        return self.duration_

    def stderr(self, max_length: int | None = None) -> str:
        payload = "".join(self.stderr_)
        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    def stdout(self, max_length: int | None = None) -> str:
        if isinstance(self.stdout_, pd.DataFrame):
            payload = str(self.stdout_)
        else:
            payload = "".join(self.stdout_)

        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    def args(self) -> List[str]:
        return self.args_

    def failures(self) -> Dict[str, str]:
        return self.failures_

    def clear(self):
        self.args_ = []
        self.returncode_ = None
        self.stdout_ = None
        self.stderr_ = None
        self.failures_ = {}
        self.title_ = None

    def status(self) -> Dict[str, Any]:
        s = {
            "title": str(self.title_),
            "args": self.args()[1:],
            "elapsed-time": str(timedelta(seconds=self.duration_)),
            "exit-code": self.returncode_,
            "notes": [],
            "status": "PASS" if len(self.failures_) == 0 else "FAIL",
        }

        if s["status"] == "PASS":
            return s

        for k, v in self.failures_.items():
            if len(v) == 0:
                s["notes"].append(k)
            else:
                s["notes"].append(f"{k}: {v}")

        return s
