# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import logging
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import is_hic


class HictkConvertCli(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-convert-cli"

    def _validate(self, expect_failure: bool):
        if expect_failure:
            if len(self.stderr()) == 0:
                self._failures["missing help message"] = ""
            if len(self.stdout()) != 0:
                self._failures["unexpected output on stdout"] = self.stdout(500).strip()
            if expect_failure and self.returncode == 0:
                self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"
            return

        if len(self.stdout()) == 0:
            self._failures["missing help message"] = ""
        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        elif not expect_failure and self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"


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

            self._failures |= validators.compare_normalizations(expected, found)

            # TODO remove once hictkpy.File().weights() is available
            # https://github.com/paulsengroup/hictkpy/pull/49
            if not is_hic(reference_file) and not is_hic(test_file):
                _, expected = self._fetch_table(reference_file, resolution=res, table="weights")
                _, found = self._fetch_table(test_file, resolution=res, table="weights")
                self._failures |= validators.compare_weights(expected, found)

    def run(  # noqa
        self,
        args: List[str],
        reference_file: str,
        test_file: str,
        resolutions: List[int],
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

        try:
            t0 = timer()
            self._run_hictk(args, timeout=timeout, env_variables=env_variables)
            t1 = timer()
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        try:
            self._validate(
                reference_file=reference_file,
                test_file=test_file,
                resolutions=resolutions,
                expect_failure=expect_failure,
            )
            t2 = timer()
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
