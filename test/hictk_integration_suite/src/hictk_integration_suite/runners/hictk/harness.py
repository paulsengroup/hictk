# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from datetime import timedelta
from timeit import default_timer as timer
from typing import Any, Dict, List, Tuple

import cooler
import pandas as pd

from hictk_integration_suite.runners import Runner
from hictk_integration_suite.validators.file_formats import is_cooler, is_hic


class HictkTestHarness:
    def __init__(
        self,
        hictk_exec: pathlib.Path | str,
        cwd: pathlib.Path | str | None = None,
        tmpdir: pathlib.Path | str | None = None,
    ):
        self._exec = pathlib.Path(hictk_exec) if hictk_exec else None
        self._cwd = pathlib.Path(cwd) if cwd else None
        self._tmpdir = pathlib.Path(tmpdir) if tmpdir else None

        self._id = None
        self._args = []
        self._expect_failure = None

        self._returncode = None
        self._stdout = None
        self._stderr = None

        self._failures = {}
        self._title = None
        self._duration = None

    def _get_hictk_keyword_option(self, option_name: str, default=None) -> Any:
        try:
            i = self._args.index(option_name)
        except ValueError:
            return default

        if i + 1 == len(self._args):
            return default

        return self._args[i + 1]

    def _get_hictk_flag_value(self, flag_name: str) -> bool:
        return flag_name in self._args

    def _run_hictk(
        self,
        args_: List[str],
        timeout: int = 1,
        env_variables: Dict[str, str] | None = None,
        colnames: List[str] | str | None = None,
    ):
        with Runner(self._exec, args_, cwd=self._cwd, tmpdir=self._tmpdir) as runner:
            self._args = runner.args
            self._returncode, self._stdout, self._stderr = runner.run(
                timeout=timeout, env_variables=env_variables, colnames=colnames
            )

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

    def _fetch_table(self, uri: str, table: str | None = None) -> Tuple[str, Any]:
        if table is None:
            table = self._get_hictk_keyword_option("--table", "pixels")

        if table == "bins":
            data = self._fetch_bins(uri)
        elif table == "chroms":
            data = self._fetch_chroms(uri)
        elif table == "pixels":
            data = self._fetch_pixels(uri)
        elif table == "normalizations":
            data = self._fetch_normalizations(uri)
        elif table == "resolutions":
            data = self._fetch_resolutions(uri)
        elif table == "cells":
            data = self._fetch_cells(uri)
        elif table == "weights":
            data = self._fetch_weights(uri)
        else:
            raise NotImplementedError

        if isinstance(data, pd.DataFrame):
            data = self._normalize_dtypes(data)

        return table, data

    def _fetch_bins(self, uri: str) -> pd.DataFrame:
        if is_hic(uri):
            return self._hictkpy_fetch_bins(uri)
        return self._cooler_fetch_bins(uri)

    def _fetch_chroms(self, uri: str) -> Dict[str, int]:
        if is_hic(uri):
            return self._hictkpy_fetch_chroms(uri)
        return self._cooler_fetch_chroms(uri)

    def _fetch_pixels(self, uri: str) -> pd.DataFrame | None:
        if is_hic(uri):
            return self._hictkpy_fetch_pixels(uri)
        return self._cooler_fetch_pixels(uri)

    def _fetch_normalizations(self, uri: str) -> List[str]:
        if is_hic(uri):
            return self._hictkpy_fetch_normalizations(uri)
        return self._cooler_fetch_normalizations(uri)

    def _fetch_resolutions(self, uri: str) -> List[int]:
        if is_hic(uri):
            return self._hictkpy_fetch_resolutions(uri)
        return self._cooler_fetch_resolutions(uri)

    def _fetch_cells(self, uri: str) -> List[str]:
        return self._cooler_fetch_cells(uri)

    def _fetch_weights(self, uri: str) -> pd.DataFrame:
        if is_hic(uri):
            return self._hictkpy_fetch_weights(uri)
        return self._cooler_fetch_weights(uri)

    def _cooler_fetch_gw_pixels(self, uri: str) -> pd.DataFrame:
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

    def _cooler_fetch_cis_only_pixels(self, uri: str) -> pd.DataFrame:
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

    def _cooler_fetch_trans_only_pixels(self, uri: str) -> pd.DataFrame:
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

    def _cooler_fetch_pixels(self, uri: str) -> pd.DataFrame | None:
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
            return self._cooler_fetch_cis_only_pixels(uri)

        if trans_only:
            # assert not cis_only
            # assert range1 is None
            # assert range2 is None
            return self._cooler_fetch_trans_only_pixels(uri)

        if range1 is None:
            # assert range2 is None
            return self._cooler_fetch_gw_pixels(uri)

        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = cooler.Cooler(uri)

        fetch_args = [arg for arg in (range1, range2) if arg is not None]
        df = clr.matrix(balance=balance, join=join, as_pixels=True).fetch(*fetch_args)
        if "balanced" in df:
            df["count"] = df["balanced"]
            df = df.drop(columns="balanced")
        return df

    def _cooler_fetch_bins(self, uri: str) -> pd.DataFrame | None:
        assert not is_hic(uri)
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

    def _cooler_fetch_chroms(self, uri: str) -> Dict[str, int]:
        assert not is_hic(uri)
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
    def _cooler_fetch_normalizations(uri: str) -> List[str]:
        assert not is_hic(uri)
        return cooler.Cooler(uri).bins().columns.drop(["chrom", "start", "end"]).tolist()

    @staticmethod
    def _cooler_fetch_resolutions(uri: str) -> List[int]:
        assert not is_hic(uri)
        path = cooler.Cooler(uri).filename
        if cooler.fileops.is_scool_file(path):
            groups = [cooler.fileops.list_coolers(path)[0]]
        elif cooler.fileops.is_multires_file(path):
            groups = cooler.fileops.list_coolers(path)
        else:
            groups = ["/"]

        return [cooler.Cooler(f"{path}::{grp}").binsize for grp in groups]

    @staticmethod
    def _cooler_fetch_cells(uri: str) -> List[str]:
        assert not is_hic(uri)
        path = cooler.Cooler(uri).filename
        if cooler.fileops.is_scool_file(path):
            return [grp.removeprefix("/cells/") for grp in cooler.fileops.list_coolers(path)]

        return []

    @staticmethod
    def _cooler_fetch_weights(uri: str) -> pd.DataFrame:
        assert not is_hic(uri)
        clr = cooler.Cooler(uri)
        columns = clr.bins().columns.drop(["chrom", "start", "end"]).tolist()
        return clr.bins()[columns][:]

    def _hictkpy_fetch_pixels(self, uri: str) -> pd.DataFrame | None:
        raise NotImplementedError

    def _hictkpy_fetch_bins(self, uri: str) -> pd.DataFrame | None:
        raise NotImplementedError

    def _hictkpy_fetch_chroms(self, uri: str) -> Dict[str, int]:
        raise NotImplementedError

    @staticmethod
    def _hictkpy_fetch_normalizations(uri: str) -> List[str]:
        raise NotImplementedError

    @staticmethod
    def _hictkpy_fetch_resolutions(uri: str) -> List[int]:
        raise NotImplementedError

    @staticmethod
    def _hictkpy_fetch_cells(uri: str) -> List[str]:
        raise NotImplementedError

    @staticmethod
    def _hictkpy_fetch_weights(uri: str) -> pd.DataFrame:
        raise NotImplementedError

    def _handle_expected_failure(self):
        if len(self.stderr()) == 0:
            self._failures["missing error message"] = ""
        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode == 0:
            self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"

    def _validate(self, expect_failure: bool):
        raise NotImplementedError

    def run(
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
        self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        self._validate(expect_failure=expect_failure)
        self._duration = timer() - t0

        return self.status()

    def ok(self) -> bool:
        return self._returncode == 0 and len(self._failures) == 0

    @property
    def returncode(self) -> int:
        return self._returncode

    @property
    def duration(self) -> float:
        return self._duration

    def stderr(self, max_length: int | None = None) -> str:
        payload = "".join(self._stderr)
        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    def stdout(self, max_length: int | None = None) -> str:
        if isinstance(self._stdout, pd.DataFrame):
            if len(self._stdout) == 0:
                payload = ""
            else:
                payload = str(self._stdout)
        else:
            payload = "".join(self._stdout)

        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[:max_length]}\n -- truncated"
        return payload

    @property
    def args(self) -> List[str]:
        return self._args

    @property
    def failures(self) -> Dict[str, str]:
        return self._failures

    def clear(self):
        self._id = None
        self._args = []
        self._expect_failure = None
        self._returncode = None
        self._stdout = None
        self._stderr = None
        self._failures = {}
        self._title = None

    def status(self) -> Dict[str, Any]:
        s = {
            "id": str(self._id),
            "title": str(self._title),
            "args": self.args[1:],
            "elapsed-time": str(timedelta(seconds=self._duration)),
            "exit-code": self._returncode,
            "expect-failure": self._expect_failure,
            "notes": [],
            "status": "PASS" if len(self._failures) == 0 else "FAIL",
        }

        if s["status"] == "PASS":
            return s

        for k, v in self._failures.items():
            if len(v) == 0:
                s["notes"].append(k)
            else:
                s["notes"].append(f"{k}: {v}")

        return s
