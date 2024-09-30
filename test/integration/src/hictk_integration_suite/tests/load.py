# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd
from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness

from .cli import HictkCli


class HictkLoadCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-load-cli"


class HictkLoad(HictkTestHarness):

    def __repr__(self) -> str:
        return "hictk-load"

    def _validate(
        self,
        test_file: pathlib.Path | str,
        reference_uri: pathlib.Path | str | None,
        resolution: int | None,
        no_validate: bool,
        expect_failure: bool,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        if no_validate:
            return

        _, found = self._fetch_table(test_file, resolution=resolution, table="bins")
        _, expected = self._fetch_table(test_file, resolution=resolution, table="bins")

        self._failures |= validators.compare_bins(found, expected)

        _, found = self._fetch_table(test_file, resolution=resolution, table="pixels")
        _, expected = self._fetch_table(reference_uri, resolution=resolution, table="pixels")
        self._failures |= validators.compare_pixels(expected, found)

    def run(  # noqa
        self,
        args: List[str],
        test_file: pathlib.Path | str,
        reference_uri: pathlib.Path | str | None,
        resolution: int | None,
        no_validate: bool,
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
            self._validate(
                test_file=test_file,
                reference_uri=reference_uri,
                resolution=resolution,
                no_validate=no_validate,
                expect_failure=expect_failure,
            )
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
