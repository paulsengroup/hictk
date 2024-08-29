# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness

from .cli import HictkCli


class HictkValidateCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-validate-cli"


class HictkValidate(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-validate"

    def _validate(self, expect_failure: bool):
        if expect_failure:
            if self.returncode != 1 and len(self.stderr()) == 0:
                self._failures["missing error message"] = ""
            if self.returncode == 1 and len(self.stdout()) == 0:
                self._failures["missing output on stdout"] = ""
            if self.returncode == 0:
                self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"
                self._failures["stdout"] = self.stdout(500)
                self._failures["stderr"] = self.stderr(500)
            return

        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stdout"] = self.stdout(500)
            self._failures["stderr"] = self.stderr(500)
            return

        output_fmt = self._get_hictk_keyword_option("--output-format", "json")
        if output_fmt not in {"json", "toml", "yaml"}:
            raise NotImplementedError

        payload = "".join(self.stdout())
        try:
            if output_fmt == "json":
                validators.metadata.json(payload)
            elif output_fmt == "toml":
                validators.metadata.toml(payload)
            elif output_fmt == "yaml":
                validators.metadata.yaml(payload)
        except (RuntimeError, ValueError) as e:
            self._failures[f"output is not valid {output_fmt}"] = str(e)

    def run(  # noqa
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
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(expect_failure=expect_failure)
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
