# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import logging
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness


class HictkDumpCli(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-metadata-dump"

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


class HictkDump(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-dump"

    def _validate(self, reference_clr: str, expect_failure: bool):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            return

        table, expected = self._fetch_table(reference_clr)
        if isinstance(self._stdout, pd.DataFrame):
            found = self._normalize_dtypes(self._stdout)
            if table == "bins":
                self._failures |= validators.compare_bins(expected, found)
            elif table == "chroms":
                self._failures |= validators.compare_chroms(expected, found)
            elif table == "pixels":
                self._failures |= validators.compare_pixels(expected, found)
            elif table == "normalizations":
                self._failures |= validators.compare_normalizations(expected, found)
            elif table == "resolutions":
                self._failures |= validators.compare_resolutions(expected, found)
            elif table == "cells":
                self._failures |= validators.compare_cells(expected, found)
            elif table == "weights":
                self._failures |= validators.compare_weights(expected, found)
            else:
                raise NotImplementedError
        else:
            self._failures["failed to read stdout into a dataframe"] = self.stdout(500).strip()

    def run(  # noqa
        self,
        args: List[str],
        reference_uri: str,
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

        table = self._get_hictk_keyword_option("--table", "pixels")

        if table == "chroms":
            colnames = ["chrom", "size"]
        elif table == "bins":
            colnames = ["chrom", "start", "end"]
        elif table == "pixels":
            if self._get_hictk_flag_value("--join"):
                colnames = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]
            else:
                colnames = ["bin1_id", "bin2_id", "count"]
        elif table == "normalizations":
            colnames = ["normalization"]
        elif table == "resolutions":
            colnames = ["resolution"]
        elif table == "cells":
            colnames = ["cell"]
        elif table == "weights":
            colnames = "infer"
        else:
            colnames = None

        t0 = timer()
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables, colnames=colnames)
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        try:
            self._validate(reference_clr=reference_uri, expect_failure=expect_failure)
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise

        self._duration = timer() - t0

        return self.status()
