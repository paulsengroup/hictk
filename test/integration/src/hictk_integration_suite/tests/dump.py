# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import os
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd
import structlog
from hictk_integration_suite import validators
from hictk_integration_suite.runners.common import normalize_df_dtypes
from hictk_integration_suite.runners.hictk import HictkTestHarness

from .cli import HictkCli


class HictkDumpCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-dump"


class HictkDump(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-dump"

    def _handle_expected_failure(self):
        if len(self.stderr()) == 0:
            self._failures["missing error message"] = ""
        if len(self.stdout()) != 0 and self.stdout() != "No columns to parse from file":
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode == 0:
            self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"

    def _validate(
        # noqa
        self,
        reference_clr: str,
        expect_failure: bool,
        table: str,
        resolution: int | None,
        range1: str | None,
        range2: str | None,
        balance: str | bool,
        cis_only: bool,
        trans_only: bool,
        join: bool,
        excluded_norms: List[str] | None,
    ):

        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            return

        _, expected = self._fetch_table(
            reference_clr,
            resolution=resolution,
            table=table,
            range1=range1,
            range2=range2,
            balance=balance,
            join=join,
            cis_only=cis_only,
            trans_only=trans_only,
        )

        if isinstance(self._stdout, pd.DataFrame):
            found = normalize_df_dtypes(self._stdout)
            if table == "bins":
                self._failures |= validators.compare_bins(expected, found)
            elif table == "chroms":
                self._failures |= validators.compare_chroms(expected, found)
            elif table == "pixels":
                self._failures |= validators.compare_pixels(expected, found)
            elif table == "normalizations":
                self._failures |= validators.compare_normalizations(expected, found, excluded_norms)
            elif table == "resolutions":
                self._failures |= validators.compare_resolutions(expected, found)
            elif table == "cells":
                self._failures |= validators.compare_cells(expected, found)
            elif table == "weights":
                self._failures |= validators.compare_weights(expected, found, excluded_norms)
            else:
                raise NotImplementedError
        elif table in {"bins", "chrom", "pixels", "resolutions"}:
            self._failures["failed to read stdout into a dataframe"] = self.stdout(500).strip()

    def run(  # noqa
        self,
        args: List[str],
        reference_uri: str,
        excluded_norms: List[str] | None = None,
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

        table = self._get_hictk_keyword_option("--table", "pixels")
        resolution = self._get_hictk_keyword_option("--resolution")
        if resolution is not None:
            resolution = int(resolution)
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")
        cis_only = self._get_hictk_flag_value("--cis-only")
        trans_only = self._get_hictk_flag_value("--trans-only")

        if table == "chroms":
            colnames = ["chrom", "size"]
        elif table == "bins":
            colnames = ["chrom", "start", "end"]
        elif table == "pixels":
            if join:
                colnames = [
                    "chrom1",
                    "start1",
                    "end1",
                    "chrom2",
                    "start2",
                    "end2",
                    "count",
                ]
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

        if env_variables is None:
            env_variables = os.environ.copy()
        else:
            env_variables = dict(env_variables.copy())

        if "LLVM_PROFILE_FILE" in env_variables:
            env_variables["LLVM_PROFILE_FILE"] = env_variables["LLVM_PROFILE_FILE"].replace("%id", str(id))

        t0 = timer()
        try:
            self._run_hictk(
                args, timeout=timeout, env_variables=env_variables, colnames=colnames, max_attempts=max_attempts
            )
        except:  # noqa
            structlog.get_logger().error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(
                reference_clr=reference_uri,
                expect_failure=expect_failure,
                table=table,
                resolution=resolution,
                range1=range1,
                range2=range2,
                balance=balance,
                join=join,
                cis_only=cis_only,
                trans_only=trans_only,
                excluded_norms=excluded_norms,
            )
        except:  # noqa
            structlog.get_logger().error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
