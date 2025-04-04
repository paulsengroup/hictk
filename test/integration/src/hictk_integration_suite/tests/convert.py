# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import os
from timeit import default_timer as timer
from typing import Any, Dict, List

import structlog
from hictk_integration_suite import validators
from hictk_integration_suite.common import URI
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import get_format, is_hic

from .cli import HictkCli


class HictkConvertCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-convert-cli"


class HictkConvert(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-convert"

    def _validate(
        self,
        test_file: str,
        reference_file: str,
        resolutions: List[int],
        expect_failure: bool,
    ):  # noqa
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

        expected_output_format = URI(test_file).path.suffix.lstrip(".")
        found_output_format = get_format(test_file)
        if expected_output_format != found_output_format:
            self._failures["unexpected output format"] = (
                f"expected {expected_output_format}, found {get_format(reference_file)}"
            )
            return

        if len(resolutions) == 0:
            expected_resolutions = self._fetch_table(reference_file, table="resolutions")[1]
        else:
            expected_resolutions = resolutions

        _, found_resolutions = self._fetch_table(test_file, table="resolutions")

        for res in expected_resolutions:
            if res not in found_resolutions:
                self._failures[f"missing data for resolution {res}"] = ""
                continue

            _, expected = self._fetch_table(reference_file, resolution=res, table="chroms")
            _, found = self._fetch_table(test_file, resolution=res, table="chroms")
            self._failures |= validators.compare_chroms(expected, found)

            _, expected = self._fetch_table(reference_file, resolution=res, table="bins")
            _, found = self._fetch_table(test_file, resolution=res, table="bins")
            self._failures |= validators.compare_bins(expected, found)

            _, expected = self._fetch_table(reference_file, resolution=res, table="pixels")
            _, found = self._fetch_table(test_file, resolution=res, table="pixels")
            self._failures |= validators.compare_pixels(expected, found)

            _, expected = self._fetch_table(reference_file, resolution=res, table="normalizations")
            _, found = self._fetch_table(test_file, resolution=res, table="normalizations")

            self._failures |= validators.compare_normalizations(expected, found, ignored_norms=["ICE", "weight"])

            if expected == found:
                _, expected = self._fetch_table(reference_file, resolution=res, table="weights")
                _, found = self._fetch_table(test_file, resolution=res, table="weights")
                self._failures |= validators.compare_weights(expected, found, ignored_weights=["ICE", "weight"])

    def run(  # noqa
        self,
        args: List[str],
        reference_file: str,
        test_file: str,
        resolutions: List[int],
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

        if env_variables is None:
            env_variables = os.environ.copy()
        else:
            env_variables = dict(env_variables.copy())

        if "LLVM_PROFILE_FILE" in env_variables:
            env_variables["LLVM_PROFILE_FILE"] = env_variables["LLVM_PROFILE_FILE"].replace("%id", str(id))

        t0 = timer()
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)
        except:  # noqa
            structlog.get_logger().error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(
                reference_file=reference_file,
                test_file=test_file,
                resolutions=resolutions,
                expect_failure=expect_failure,
            )
        except:  # noqa
            structlog.get_logger().error(f"failed to validate output produced by {args}")
            raise

        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
