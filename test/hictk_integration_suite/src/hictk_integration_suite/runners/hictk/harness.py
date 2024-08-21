# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from datetime import timedelta
from timeit import default_timer as timer
from typing import Any, Dict, List, Tuple

import cooler
import hictkpy
import pandas as pd

from hictk_integration_suite.runners import Runner
from hictk_integration_suite.validators.file_formats import (
    is_cooler,
    is_hic,
    is_multires,
)


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

    @staticmethod
    def _filter_bins(df: pd.DataFrame, chroms: Dict[str, int], range1: str | None, range2: str | None) -> pd.DataFrame:
        if range1 is None:
            # assert range2 is None
            return df

        chrom1, start1, end1 = cooler.api.parse_region(range1, chroms)
        df1 = df[(df["chrom"] == chrom1) & (df["start"] >= start1) & (df["start"] < end1)]
        if range2 is None or range1 == range2:
            return df1

        chrom2, start2, end2 = cooler.api.parse_region(range2, chroms)
        df2 = df[(df["chrom"] == chrom2) & (df["start"] >= start2) & (df["start"] < end2)]

        return pd.concat([df1, df2])

    @staticmethod
    def _filter_chroms(chroms: Dict[str, int], range1: str | None, range2: str | None) -> Dict[str, int]:
        if range1 is None:
            # assert range2 is None
            return chroms

        chrom1, _, _ = range1.partition(":")
        if range2 is None:
            return {chrom1: chroms[chrom1]}

        chrom2, _, _ = range2.partition(":")
        return {chrom1: chroms[chrom1], chrom2: chroms[chrom2]}

    @staticmethod
    def _filter_weights(
        df: pd.DataFrame, chroms: Dict[str, int], range1: str | None, range2: str | None
    ) -> pd.DataFrame:
        columns = df.columns.drop(["chrom", "start", "end"]).tolist()
        if range1 is None:
            # assert range2 is None
            return df[columns]

        chrom1, start1, end1 = cooler.api.parse_region(range1, chroms)
        df1 = df[(df["chrom"] == chrom1) & (df["start"] >= start1) & (df["start"] < end1)]
        if range2 is None or range1 == range2:
            return df1[columns]

        chrom2, start2, end2 = cooler.api.parse_region(range2, chroms)
        df2 = df[(df["chrom"] == chrom2) & (df["start"] >= start2) & (df["start"] < end2)]

        return pd.concat([df1, df2])[columns]

    def _fetch_table(self, uri: str, table: str | None = None) -> Tuple[str, Any]:
        if table is None:
            table = self._get_hictk_keyword_option("--table", "pixels")

        resolution = self._get_hictk_keyword_option("--resolution")
        path, _, _ = uri.partition("::")

        if table == "bins":
            data = self._fetch_bins(path, resolution)
        elif table == "chroms":
            data = self._fetch_chroms(path)
        elif table == "pixels":
            data = self._fetch_pixels(path, resolution)
        elif table == "normalizations":
            data = self._fetch_normalizations(path, resolution)
        elif table == "resolutions":
            data = self._fetch_resolutions(path)
        elif table == "cells":
            data = self._fetch_cells(path)
        elif table == "weights":
            data = self._fetch_weights(path, resolution)
        else:
            raise NotImplementedError

        if isinstance(data, pd.DataFrame):
            data = self._normalize_dtypes(data)

        return table, data

    def _fetch_bins(self, path: str, resolution: int | None) -> pd.DataFrame:
        if is_hic(path):
            return self._hictkpy_fetch_bins(path, resolution)
        return self._cooler_fetch_bins(path, resolution)

    def _fetch_chroms(self, path: str) -> Dict[str, int]:
        if is_hic(path):
            return self._hictkpy_fetch_chroms(path)
        return self._cooler_fetch_chroms(path)

    def _fetch_pixels(self, path: str, resolution: int | None) -> pd.DataFrame | None:
        if is_hic(path):
            return self._hictkpy_fetch_pixels(path, resolution)
        return self._cooler_fetch_pixels(path, resolution)

    def _fetch_normalizations(self, path: str, resolution: int | None) -> List[str]:
        if is_hic(path):
            return self._hictkpy_fetch_normalizations(path, resolution)
        return self._cooler_fetch_normalizations(path, resolution)

    def _fetch_resolutions(self, path: str) -> List[int]:
        if is_hic(path):
            return self._hictkpy_fetch_resolutions(path)
        return self._cooler_fetch_resolutions(path)

    def _fetch_cells(self, path: str) -> List[str]:
        return self._cooler_fetch_cells(path)

    def _fetch_weights(self, path: str, resolution: int | None) -> pd.DataFrame:
        if is_hic(path):
            return self._hictkpy_fetch_weights(path, resolution)
        return self._cooler_fetch_weights(path, resolution)

    @staticmethod
    def _open_cooler(uri: str, resolution: int | None) -> cooler.Cooler | None:
        try:
            if resolution is None:
                return cooler.Cooler(uri)
            if is_cooler(uri):
                clr = cooler.Cooler(uri)
                assert clr.binsize == resolution
                return clr
            return cooler.Cooler(f"{uri}::/resolutions{resolution}")
        except:  # noqa
            return None

    def _cooler_fetch_gw_pixels(self, path: str, resolution: int | None) -> pd.DataFrame | None:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = self._open_cooler(path, resolution)
        if not clr:
            return None
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

    def _cooler_fetch_cis_only_pixels(self, path: str, resolution: int | None) -> pd.DataFrame | None:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = self._open_cooler(path, resolution)
        if not clr:
            return None
        sel = clr.matrix(balance=balance, join=join, as_pixels=True)

        dfs = []
        for i1, chrom1 in enumerate(clr.chromnames):
            df = sel.fetch(chrom1)
            if "balanced" in df:
                df["count"] = df["balanced"]
                df = df.drop(columns="balanced")
            dfs.append(df)

        return pd.concat(dfs)

    def _cooler_fetch_trans_only_pixels(self, path: str, resolution: int | None) -> pd.DataFrame:
        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = self._open_cooler(path, resolution)
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

    def _cooler_fetch_pixels(self, path: str, resolution: int | None = None) -> pd.DataFrame | None:
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")
        cis_only = self._get_hictk_flag_value("--cis-only")
        trans_only = self._get_hictk_flag_value("--trans-only")

        if cis_only:
            # assert not trans_only
            # assert range1 is None
            # assert range2 is None
            return self._cooler_fetch_cis_only_pixels(path, resolution)

        if trans_only:
            # assert not cis_only
            # assert range1 is None
            # assert range2 is None
            return self._cooler_fetch_trans_only_pixels(path, resolution)

        if range1 is None:
            # assert range2 is None
            return self._cooler_fetch_gw_pixels(path, resolution)

        balance = self._get_hictk_keyword_option("--balance", False)
        join = self._get_hictk_flag_value("--join")

        clr = self._open_cooler(path, resolution)
        if not clr:
            return None

        fetch_args = [arg for arg in (range1, range2) if arg is not None]
        df = clr.matrix(balance=balance, join=join, as_pixels=True).fetch(*fetch_args)
        if "balanced" in df:
            df["count"] = df["balanced"]
            df = df.drop(columns="balanced")
        return df

    def _cooler_fetch_bins(self, path: str, resolution: int | None = None) -> pd.DataFrame | None:
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        clr = self._open_cooler(path, resolution)
        return self._filter_bins(clr.bins()[["chrom", "start", "end"]][:], clr.chromsizes.to_dict(), range1, range2)

    def _cooler_fetch_chroms(self, path: str) -> Dict[str, int]:
        uri = path
        if not cooler.fileops.is_cooler(path):
            # Try to open a .scool or .mcool
            groups = cooler.fileops.list_coolers(path)
            if len(groups) == 0:
                raise RuntimeError(f'unable to find any cooler files under "{path}"')

            uri = f"{path}::{groups[0]}"

        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        clr = self._open_cooler(uri, None)
        if not clr:
            raise RuntimeError(f'unable to fetch chromosomes from "{path}"')
        return self._filter_chroms(clr.chromsizes.to_dict(), range1, range2)

    def _cooler_fetch_normalizations(self, path: str, resolution: int | None = None) -> List[str]:
        clr = self._open_cooler(path, resolution)
        if not clr:
            raise RuntimeError(f'unable to fetch normalizations from "{path}"')
        return clr.bins().columns.drop(["chrom", "start", "end"]).tolist()

    def _cooler_fetch_resolutions(self, path: str, resolution: int | None) -> List[int]:
        clr = self._open_cooler(path, resolution)
        if clr:
            return [clr.binsize]

        if cooler.fileops.is_scool_file(path):
            groups = [cooler.fileops.list_coolers(path)[0]]
        elif cooler.fileops.is_multires_file(path):
            groups = cooler.fileops.list_coolers(path)
        else:
            raise RuntimeError(f'unable to fetch resolutions from "{path}"')

        return [cooler.Cooler(f"{path}::{grp}").binsize for grp in groups]

    @staticmethod
    def _cooler_fetch_cells(path: str) -> List[str]:
        if cooler.fileops.is_scool_file(path):
            return [grp.removeprefix("/cells/") for grp in cooler.fileops.list_coolers(path)]
        return []

    def _cooler_fetch_weights(self, path: str, resolution: int | None) -> pd.DataFrame:
        clr = self._open_cooler(path, resolution)
        if not clr:
            raise RuntimeError(f'unable to fetch weights from "{path}"')

        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        return self._filter_weights(clr.bins()[:], clr.chromsizes.to_dict(), range1, range2)

    @staticmethod
    def _hictkpy_fetch_cis_only_pixels(
        path: str, resolution: int | None, normalization: str, join: bool
    ) -> pd.DataFrame:
        f = hictkpy.File(path, resolution)

        chromnames = list(f.chroms().keys())
        dfs = []
        for i1, chrom1 in enumerate(chromnames):
            dfs.append(f.fetch(chrom1, normalization=normalization, join=join).to_df())

        return pd.concat(dfs)

    @staticmethod
    def _hictkpy_fetch_trans_only_pixels(
        path: str, resolution: int | None, normalization: str, join: bool
    ) -> pd.DataFrame:
        f = hictkpy.File(path, resolution)

        # TODO: optimize
        df = f.fetch(normalization=normalization, join=join).to_df()
        return df[df["chrom1"] != df["chrom2"]].reset_index(drop=True)

    def _hictkpy_fetch_pixels(self, path, resolution: int | None = None) -> pd.DataFrame | None:
        if resolution is None:
            resolution = self._get_hictk_keyword_option("--resolution")

        range1 = self._get_hictk_keyword_option("--range", "")
        range2 = self._get_hictk_keyword_option("--range2", "")
        cis_only = self._get_hictk_flag_value("--cis-only")
        trans_only = self._get_hictk_flag_value("--trans-only")
        normalization = self._get_hictk_keyword_option("--balance", "NONE")
        join = self._get_hictk_flag_value("--join")

        if cis_only:
            # assert not trans_only
            # assert range1 is None
            # assert range2 is None
            return self._hictkpy_fetch_cis_only_pixels(path, resolution, normalization, join)

        if trans_only:
            # assert not cis_only
            # assert range1 is None
            # assert range2 is None
            return self._hictkpy_fetch_trans_only_pixels(path, resolution, normalization, join)

        return hictkpy.File(path, resolution).fetch(range1, range2, normalization=normalization, join=join).to_df()

    def _hictkpy_fetch_bins(self, path: str, resolution: int | None) -> pd.DataFrame | None:
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        return self._filter_bins(hictkpy.File(path, resolution).bins(), range1, range2)

    def _hictkpy_fetch_chroms(self, path: str) -> Dict[str, int]:
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")
        if is_multires(path):
            f = hictkpy.MultiResFile(path)
        else:
            f = hictkpy.File(path)
        return self._filter_chroms(f.chroms(), range1, range2)

    @staticmethod
    def _hictkpy_fetch_normalizations(path: str) -> List[str]:
        resolution = None
        if is_multires(path):
            resolution = hictkpy.MultiResFile(path).resolutions()[-1]
        return hictkpy.File(path, resolution).avail_normalizations()

    @staticmethod
    def _hictkpy_fetch_resolutions(path) -> List[int]:
        if is_multires(path):
            return hictkpy.MultiResFile(path).resolutions()
        return [hictkpy.File(path).resolution()]

    def _hictkpy_fetch_weights(self, path: str, resolution: int | None) -> pd.DataFrame:
        range1 = self._get_hictk_keyword_option("--range")
        range2 = self._get_hictk_keyword_option("--range2")

        f = hictkpy.File(path, resolution)
        data = {}
        for norm in f.avail_normalizations():
            data[norm] = f.weights(norm)

        return self._filter_weights(pd.DataFrame(data), f.chroms(), range1, range2)

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
