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


class HictkMergeCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-merge-cli"


class HictkMerge(HictkTestHarness):

    def __repr__(self) -> str:
        return "hictk-merge"

    def _merge_pixels(self, dfs: List[pd.DataFrame]) -> pd.DataFrame:
        if len(dfs) == 1:
            return dfs[0]

        assert "bin1_id" in dfs[0].columns
        assert "bin2_id" in dfs[0].columns
        assert "count" in dfs[0].columns

        df = pd.concat([dff[["bin1_id", "bin2_id"]] for dff in dfs]).drop_duplicates(keep="first")
        if pd.api.types.is_integer_dtype(dfs[0]["count"].dtype):
            df["count"] = 0
        else:
            df["count"] = 0.0

        for dff in dfs:
            df = df.merge(dff, on=["bin1_id", "bin2_id"], how="left", suffixes=("", "_new"))
            df["count"] += df["count_new"]
            df = df.drop(columns="count_new")

        return df

    def _validate(
        self,
        input_files: List[pathlib.Path | str],
        test_file: pathlib.Path | str,
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

        resolution = self._get_hictk_keyword_option("--resolution")
        if resolution is not None:
            resolution = int(resolution)
        expected = self._merge_pixels([self._fetch_table(f, resolution, table="pixels")[1] for f in input_files])
        found = self._fetch_table(test_file, resolution, table="pixels")[1]

        self._failures = validators.compare_pixels(expected, found)

    def run(  # noqa
        self,
        args: List[str],
        input_files: List[pathlib.Path | str],
        test_file: pathlib.Path | str,
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
            self._validate(
                input_files=input_files,
                test_file=test_file,
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
