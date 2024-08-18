# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from datetime import timedelta
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite.runners import Runner


class HictkTestHarness:
    def __init__(self, hictk_exec: pathlib.Path, cwd: str | None = None, tmpdir: pathlib.Path | None = None):
        self.exec = hictk_exec
        self.cwd = cwd
        self.tmpdir = tmpdir

        self.args_ = []

        self.returncode_ = None
        self.stdout_ = None
        self.stderr_ = None

        self.failures_ = {}
        self.title_ = None
        self.duration_ = None

    def _run_hictk(self, args: List[str], timeout: int = 3600, env_variables: Dict[str, str] | None = None):
        runner = Runner(self.exec, args)
        self.args_ = runner.args()
        self.returncode_, self.stdout_, self.stderr_ = runner.run(timeout=timeout, env_variables=env_variables)

    def _validate(self, expect_failure: bool):
        raise NotImplemented()

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
            "args": self.args(),
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
