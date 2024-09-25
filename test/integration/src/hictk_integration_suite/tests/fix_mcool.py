# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import json
import logging
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite.runners import Runner
from hictk_integration_suite.runners.hictk import HictkTestHarness

from .cli import HictkCli


class HictkFixMcoolCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-fix-mcool-cli"


class HictkFixMcool(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-fix-mcool"

    def _validate(self, test_file: pathlib.Path | str, expect_failure: bool):
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = (
                f"expected zero, found {self.returncode}; excerpt from stderr: {self.stderr(500).strip()}"
            )
            return

        return_code_, stdout_, stderr_ = Runner(
            self._exec, ["validate", str(test_file), "--validate-index", "--exhaustive"], self._cwd, self._tmpdir
        ).run(timeout=180.0)
        if return_code_ != 0 and len(stdout_) == 0:
            self._failures["mcool file is still corrupted"] = f'hictk validate failed unexpectedly: "{stderr_}"'
            return

        status = json.loads("".join(stdout_))
        if not status["is_valid_mcool"]:
            self._failures["mcool file is still corrupted"] = ""

    def run(  # noqa
        self,
        args: List[str],
        test_file: pathlib.Path | str,
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        max_attempts: int = 1,
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
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(test_file=test_file, expect_failure=expect_failure)
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
