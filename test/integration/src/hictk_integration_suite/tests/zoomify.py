# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite import validators
from hictk_integration_suite.runners.cooler import CoolerCoarsen
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import get_format, is_cooler

from .cli import HictkCli


class HictkZoomifyCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-zoomify-cli"


class HictkZoomify(HictkTestHarness):

    def __repr__(self) -> str:
        return "hictk-zoomify"

    def _validate(
        self,
        test_file: str,
        reference_file: str,
        expect_failure: bool,
        resolutions: List[int],
        expect_single_resolution: bool,
        tmpdir: pathlib.Path,
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

        assert len(resolutions) != 0

        if expect_single_resolution:
            if test_file.endswith("cool") and not is_cooler(test_file):
                self._failures["unexpected file format"] = (
                    f"expected file in .cool format, found {get_format(test_file)}"
                )

            num_avail_res = len(self._fetch_table(test_file, table="resolutions")[1])
            if num_avail_res != 1:
                self._failures["file is not single-resolution"] = (
                    f"expected to find exactly one resolution, found {num_avail_res}"
                )

        for res in resolutions:
            _, found = self._fetch_table(test_file, resolution=res, table="bins")

            if found is None:
                self._failures[f"missing resolution {res}"] = ""
                continue

            reference_uri = reference_file
            try:
                _, expected = self._fetch_table(reference_uri, resolution=res, table="bins")
            except AssertionError:
                expected = None

            if expected is None:
                logging.debug(f"coarsening cooler at URI {reference_uri} to resolution {res}...")
                reference_uri = CoolerCoarsen(
                    input_uri=reference_uri,
                    base_resolution=None,
                    output_uri=tmpdir / f"tmp.{res}.cool",
                    target_resolution=res,
                ).coarsen()
                _, expected = self._fetch_table(reference_uri, resolution=res, table="bins")
                assert expected is not None

            self._failures |= validators.compare_bins(found, expected)

            _, found = self._fetch_table(test_file, resolution=res, table="pixels")
            _, expected = self._fetch_table(reference_uri, resolution=res, table="pixels")
            self._failures |= validators.compare_pixels(expected, found)

    def run(  # noqa
        self,
        args: List[str],
        reference_file: pathlib.Path | str,
        test_file: pathlib.Path | str,
        resolutions: List[int],
        expect_single_resolution: True,
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
            tmpdir = self._get_hictk_keyword_option("--tmpdir")
            assert tmpdir is not None
            self._validate(
                reference_file=reference_file,
                test_file=test_file,
                expect_failure=expect_failure,
                resolutions=resolutions,
                expect_single_resolution=expect_single_resolution,
                tmpdir=pathlib.Path(tmpdir),
            )
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
