# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import sys
from typing import Dict, List, Tuple

import cooler
import pandas as pd

from hictk_integration_suite.runners.common import (
    filter_bins,
    filter_chroms,
    filter_weights,
    normalize_df_dtypes,
)
from hictk_integration_suite.validators.file_formats import (
    is_cooler,
    is_multires,
    is_scool,
)


class CoolerDump:
    def __init__(
        self,
        uri: pathlib.Path,
        resolution: int | None,
    ):
        self._file = pathlib.Path(str(uri).partition("::")[0])
        self._uri = str(uri)
        self._clr = self._open(str(uri), resolution)

    @staticmethod
    def _open(uri: str, resolution: int | None) -> cooler.Cooler | None:
        if is_cooler(uri):
            clr = cooler.Cooler(uri)
            if resolution is not None:
                assert clr.binsize == resolution
            return clr
        if is_multires(uri) and resolution is not None:
            return cooler.Cooler(f"{uri}::/resolutions{resolution}")
        return None

    def _fetch_gw_pixels(self, balance: str | bool, join: bool) -> pd.DataFrame | None:
        if self._clr is None:
            return None
        sel = self._clr.matrix(balance=balance, join=join, as_pixels=True, ignore_index=False)

        dfs = []
        for i1, chrom1 in enumerate(self._clr.chromnames):
            for chrom2 in self._clr.chromnames[i1:]:
                df = sel.fetch(chrom1, chrom2)
                if "balanced" in df:
                    df["count"] = df["balanced"]
                    df = df.drop(columns="balanced")
                dfs.append(df)

        df = pd.concat(dfs)
        return df.sort_index().reset_index(drop=True)

    def _fetch_cis_only_pixels(self, balance: str | bool, join: bool) -> pd.DataFrame | None:
        if self._clr is None:
            return None
        sel = self._clr.matrix(balance=balance, join=join, as_pixels=True)

        dfs = []
        for chrom in self._clr.chromnames:
            df = sel.fetch(chrom)
            if "balanced" in df:
                df["count"] = df["balanced"]
                df = df.drop(columns="balanced")

            dfs.append(df)

        return pd.concat(dfs)

    def _fetch_trans_only_pixels(self, balance: str | bool, join: bool) -> pd.DataFrame | None:
        if self._clr is None:
            return None

        sel = self._clr.matrix(balance=balance, join=join, as_pixels=True, ignore_index=False)

        dfs = []
        for i1, chrom1 in enumerate(self._clr.chromnames):
            for chrom2 in self._clr.chromnames[i1 + 1 :]:
                df = sel.fetch(chrom1, chrom2)
                if "balanced" in df:
                    df["count"] = df["balanced"]
                    df = df.drop(columns="balanced")
                dfs.append(df)

        return pd.concat(dfs).sort_index().reset_index(drop=True)

    def _fetch_pixels(
        self,
        range1: str | None,
        range2: str | None,
        balance: str | bool,
        join: bool,
        cis_only: bool,
        trans_only: bool,
    ) -> pd.DataFrame | None:
        if cis_only:
            assert not trans_only
            assert range1 is None
            assert range2 is None
            return self._fetch_cis_only_pixels(balance, join)

        if trans_only:
            assert not cis_only
            assert range1 is None
            assert range2 is None
            return self._fetch_trans_only_pixels(balance, join)

        if range1 is None:
            assert range2 is None
            return self._fetch_gw_pixels(balance, join)

        if self._clr is None:
            return None
        fetch_args = [arg for arg in (range1, range2) if arg is not None]
        df = self._clr.matrix(balance=balance, join=join, as_pixels=True).fetch(*fetch_args)
        if "balanced" in df:
            df["count"] = df["balanced"]
            df = df.drop(columns="balanced")
        return df

    def _fetch_bins(self, range1: str | None, range2: str | None) -> pd.DataFrame | None:
        clr = self._clr
        if clr is None and is_scool(self._uri):
            groups = cooler.fileops.list_coolers(str(self._file))
            if len(groups) == 0:
                raise RuntimeError(f'unable to find any cooler files under "{self._file}"')

            uri = f"{self._file}::{groups[0]}"
            clr = self._open(uri, None)
            if clr is None:
                return None

        return filter_bins(clr.bins()[["chrom", "start", "end"]][:], clr.chromsizes.to_dict(), range1, range2)

    def _fetch_chroms(self, range1: str | None, range2: str | None) -> Dict[str, int] | None:
        if not cooler.fileops.is_cooler(str(self._file)):
            # Try to open a .scool or .mcool
            groups = cooler.fileops.list_coolers(str(self._file))
            if len(groups) == 0:
                raise RuntimeError(f'unable to find any cooler files under "{self._file}"')

            uri = f"{self._file}::{groups[0]}"
            return filter_chroms(self._open(uri, None).chromsizes.to_dict(), range1, range2)

        if not self._clr:
            return None
        return filter_chroms(self._clr.chromsizes.to_dict(), range1, range2)

    def _fetch_normalizations(self) -> List[str] | None:
        if cooler.fileops.is_multires_file(self._uri):
            uri = f"{self._uri}::{cooler.fileops.list_coolers(self._uri)[-1]}"
            return cooler.Cooler(uri).bins().columns.drop(["chrom", "start", "end"]).tolist()

        if self._clr is None:
            return None
        return self._clr.bins().columns.drop(["chrom", "start", "end"]).tolist()

    def _fetch_resolutions(self) -> List[int] | None:
        if self._clr is not None:
            return [self._clr.binsize]

        if cooler.fileops.is_scool_file(self._uri):
            groups = [cooler.fileops.list_coolers(self._uri)[0]]
        elif cooler.fileops.is_multires_file(self._uri):
            groups = cooler.fileops.list_coolers(self._uri)
        else:
            return None

        return [cooler.Cooler(f"{self._uri}::{grp}").binsize for grp in groups]

    def _fetch_cells(self) -> List[str] | None:
        if cooler.fileops.is_scool_file(self._uri):
            return [grp.removeprefix("/cells/") for grp in cooler.fileops.list_coolers(self._uri)]
        return None

    def _fetch_weights(self, range1: str | None, range2: str | None) -> pd.DataFrame | None:
        if not self._clr:
            return None
        return filter_weights(self._clr.bins()[:], self._clr.chromsizes.to_dict(), range1, range2)

    def dump(
        self,
        **kwargs,
    ) -> Tuple[str, pd.DataFrame | List[int] | List[str] | None]:
        table = kwargs.get("table", "pixels")
        range1 = kwargs.get("range1")
        range2 = kwargs.get("range2")
        balance = kwargs.get("balance", False)
        join = kwargs.get("join", False)
        cis_only = kwargs.get("cis_only", False)
        trans_only = kwargs.get("trans_only", False)

        if table == "bins":
            data = self._fetch_bins(range1, range2)
        elif table == "chroms":
            data = self._fetch_chroms(range1, range2)
        elif table == "pixels":
            data = self._fetch_pixels(range1, range2, balance, join, cis_only, trans_only)
        elif table == "normalizations":
            data = self._fetch_normalizations()
        elif table == "resolutions":
            data = self._fetch_resolutions()
        elif table == "cells":
            data = self._fetch_cells()
        elif table == "weights":
            data = self._fetch_weights(range1, range2)
        else:
            raise NotImplementedError

        if isinstance(data, pd.DataFrame):
            data = normalize_df_dtypes(data)

        return table, data
