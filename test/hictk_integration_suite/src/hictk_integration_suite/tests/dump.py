# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import logging
from timeit import default_timer as timer
from typing import Any, Dict, List, Tuple

import cooler
import pandas as pd

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import is_cooler, is_hic


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

    @staticmethod
    def _normalize_dtypes(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        if "chrom1" in df.columns:
            df[["chrom1", "chrom2"]] = df[["chrom1", "chrom2"]].astype(str)
        if "chrom" in df.columns:
            df["chrom"] = df["chrom"].astype(str)

        columns = df.select_dtypes(include=int).columns.tolist()
        df[columns] = df[columns].astype(int)

        columns = df.select_dtypes(include=float).columns.tolist()
        df[columns] = df[columns].astype(float)
        return df

    def _cooler_dump(self, uri: str) -> Tuple[str, Any]:
        assert not is_hic(uri)
        table = self._get_hictk_keyword_option("--table")
        if table is None:
            table = "pixels"

        if table == "bins":
            data = self._fetch_bins_cooler(uri)
        elif table == "chroms":
            data = self._fetch_chroms_cooler(uri)
        elif table == "pixels":
            data = self._run_query_cooler(uri)
        elif table == "normalizations":
            data = self._fetch_normalizations_cooler(uri)
        elif table == "resolutions":
            data = self._fetch_resolutions_cooler(uri)
        elif table == "cells":
            data = self._fetch_cells_cooler(uri)
        elif table == "weights":
            data = self._fetch_weights_cooler(uri)
        else:
            raise NotImplementedError

        if isinstance(data, pd.DataFrame):
            data = self._normalize_dtypes(data)

        return table, data

    def _run_gw_query_cooler(self, uri: str) -> pd.DataFrame:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = cooler.Cooler(uri)
        sel = clr.matrix(balance=balance, join=join, as_pixels=True, ignore_index=False)

        dfs = []
        for i1, chrom1 in enumerate(clr.chromnames):
            for chrom2 in clr.chromnames[i1:]:
                df = sel.fetch(chrom1, chrom2)
                if "balanced" in df:
                    df["count"] = df["balanced"]
                    df = df.drop(columns="balanced")
                dfs.append(df)

        df = pd.concat(dfs)
        return df.sort_index().reset_index(drop=True)

    def _run_cis_only_query_cooler(self, uri: str) -> pd.DataFrame:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = cooler.Cooler(uri)
        sel = clr.matrix(balance=balance, join=join, as_pixels=True)

        dfs = []
        for i1, chrom1 in enumerate(clr.chromnames):
            df = sel.fetch(chrom1)
            if "balanced" in df:
                df["count"] = df["balanced"]
                df = df.drop(columns="balanced")
            dfs.append(df)

        return pd.concat(dfs)

    def _run_trans_only_query_cooler(self, uri: str) -> pd.DataFrame:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = cooler.Cooler(uri)
        sel = clr.matrix(balance=balance, join=join, as_pixels=True, ignore_index=False)

        dfs = []
        for i1, chrom1 in enumerate(clr.chromnames):
            for chrom2 in clr.chromnames[i1 + 1 :]:
                df = sel.fetch(chrom1, chrom2)
                if "balanced" in df:
                    df["count"] = df["balanced"]
                    df = df.drop(columns="balanced")
                dfs.append(df)

        return pd.concat(dfs).sort_index().reset_index(drop=True)

    def _run_query_cooler(self, uri: str) -> pd.DataFrame | None:
        if not is_cooler(uri):
            resolution = self._get_hictk_keyword_option("--resolution")
            uri = f"{uri}::/resolutions/{resolution}"
            if not is_cooler(uri):
                return None

        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")
        cis_only = self._get_hictk_flag_value("--cis-only")
        trans_only = self._get_hictk_flag_value("--trans-only")

        if cis_only:
            # assert not trans_only
            # assert range1 is None
            # assert range2 is None
            return self._run_cis_only_query_cooler(uri)

        if trans_only:
            # assert not cis_only
            # assert range1 is None
            # assert range2 is None
            return self._run_trans_only_query_cooler(uri)

        if range1 is None:
            # assert range2 is None
            return self._run_gw_query_cooler(uri)

        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = cooler.Cooler(uri)

        fetch_args = [arg for arg in (range1, range2) if arg is not None]
        df = clr.matrix(balance=balance, join=join, as_pixels=True).fetch(*fetch_args)
        if "balanced" in df:
            df["count"] = df["balanced"]
            df = df.drop(columns="balanced")
        return df

    def _fetch_bins_cooler(self, uri: str) -> pd.DataFrame | None:
        if cooler.fileops.is_multires_file(uri):
            return None
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        clr = cooler.Cooler(uri)
        df = clr.bins()[:][["chrom", "start", "end"]]
        if range1 is None:
            # assert range2 is None
            return df

        chrom1, start1, end1 = cooler.api.parse_region(range1, clr.chromsizes)
        df1 = df[(df["chrom"] == chrom1) & (df["start"] >= start1) & (df["start"] < end1)]
        if range2 is None or range1 == range2:
            return df1

        chrom2, start2, end2 = cooler.api.parse_region(range2, clr.chromsizes)
        df2 = df[(df["chrom"] == chrom2) & (df["start"] >= start2) & (df["start"] < end2)]

        return pd.concat([df1, df2])

    def _fetch_chroms_cooler(self, uri: str) -> Dict[str, int]:
        if not cooler.fileops.is_cooler(uri):
            # Try to open a .scool or .mcool
            groups = cooler.fileops.list_coolers(uri)
            if len(groups) == 0:
                raise RuntimeError(f'unable to find any cooler files under "{uri}"')

            uri = f"{uri}::{groups[0]}"

        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        clr = cooler.Cooler(uri)
        chroms = clr.chromsizes.to_dict()

        if range1 is None:
            # assert range2 is None
            return chroms

        chrom1, _, _ = cooler.api.parse_region(range1, clr.chromsizes)
        if range2 is None:
            return {chrom1: chroms[chrom1]}

        chrom2, _, _ = cooler.api.parse_region(range2, clr.chromsizes)
        return {chrom1: chroms[chrom1], chrom2: chroms[chrom2]}

    @staticmethod
    def _fetch_normalizations_cooler(uri: str) -> List[str]:
        return cooler.Cooler(uri).bins().columns.drop(["chrom", "start", "end"]).tolist()

    @staticmethod
    def _fetch_resolutions_cooler(uri: str) -> List[int]:
        path = cooler.Cooler(uri).filename
        if cooler.fileops.is_scool_file(path):
            groups = [cooler.fileops.list_coolers(path)[0]]
        elif cooler.fileops.is_multires_file(path):
            groups = cooler.fileops.list_coolers(path)
        else:
            groups = ["/"]

        return [cooler.Cooler(f"{path}::{grp}").binsize for grp in groups]

    @staticmethod
    def _fetch_cells_cooler(uri: str) -> List[str]:
        path = cooler.Cooler(uri).filename
        if cooler.fileops.is_scool_file(path):
            return [grp.removeprefix("/cells/") for grp in cooler.fileops.list_coolers(path)]

        return []

    @staticmethod
    def _fetch_weights_cooler(uri: str) -> pd.DataFrame:
        clr = cooler.Cooler(uri)
        columns = clr.bins().columns.drop(["chrom", "start", "end"]).tolist()
        return clr.bins()[columns][:]

    def _validate(self, reference_clr: str, expect_failure: bool):  # noqa
        if expect_failure:
            if len(self.stderr()) == 0:
                self._failures["missing error message"] = ""
            if len(self.stdout()) != 0:
                self._failures["unexpected output on stdout"] = self.stdout(100).strip()
            if self.returncode == 0:
                self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"
            return

        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            return

        table, expected = self._cooler_dump(reference_clr)
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
