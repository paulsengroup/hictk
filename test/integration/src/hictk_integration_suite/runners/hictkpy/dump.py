# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Dict, List, Tuple

import hictkpy
import pandas as pd

from hictk_integration_suite.runners.common import (
    filter_bins,
    filter_chroms,
    filter_weights,
    normalize_df_dtypes,
)
from hictk_integration_suite.validators.file_formats import (
    is_cooler,
    is_hic,
    is_multires,
    is_scool,
)


class HictkpyDump:
    def __init__(
        self,
        uri: pathlib.Path,
        resolution: int | None,
    ):
        self._f = self._open(str(uri), resolution)

    def __bool__(self) -> bool:
        return self._f is not None

    @staticmethod
    def _open(
        uri: str, resolution: int | None
    ) -> hictkpy.File | hictkpy.MultiResFile | hictkpy.cooler.SingleCellFile | None:
        if is_multires(uri) and resolution is None:
            return hictkpy.MultiResFile(uri)

        if is_scool(uri):
            f = hictkpy.cooler.SingleCellFile(uri)
            if resolution is not None:
                assert f.resolution() == resolution
            return f

        try:
            f = hictkpy.File(uri, resolution)
            if resolution is not None:
                assert f.resolution() == resolution
            return f
        except RuntimeError:
            return None

    def _is_single_res_file(self) -> bool:
        return isinstance(self._f, hictkpy.File)

    def _is_multi_res_file(self) -> bool:
        return isinstance(self._f, hictkpy.MultiResFile)

    def _is_single_cell_file(self) -> bool:
        return isinstance(self._f, hictkpy.cooler.SingleCellFile)

    def _fetch_cis_only_pixels(self, normalization: str, join: bool) -> pd.DataFrame | None:
        if self._f is None:
            return None

        dfs = []
        for chrom in self._f.chromosomes().keys():
            dfs.append(self._f.fetch(chrom, normalization=normalization, join=join).to_df())

        return pd.concat(dfs)

    def _fetch_trans_only_pixels(self, normalization: str, join: bool) -> pd.DataFrame | None:
        if self._f is None:
            return None

        # TODO: optimize
        df = self._f.fetch(normalization=normalization, join=join).to_df()
        return df[df["chrom1"] != df["chrom2"]].reset_index(drop=True)

    def _fetch_pixels(
        self,
        range1: str | None,
        range2: str | None,
        balance: str | bool,
        join: bool,
        cis_only: bool,
        trans_only: bool,
    ) -> pd.DataFrame | None:
        if isinstance(balance, bool):
            normalization = "weight" if balance else "NONE"
        else:
            normalization = balance

        if cis_only:
            assert not trans_only
            assert range1 is None
            assert range2 is None
            return self._fetch_cis_only_pixels(normalization, join)

        if trans_only:
            assert not cis_only
            assert range1 is None
            assert range2 is None
            return self._fetch_trans_only_pixels(normalization, join)

        if self._f is None:
            return None

        if range1 is None:
            range1 = ""
        if range2 is None:
            range2 = ""

        return self._f.fetch(range1, range2, normalization=normalization, join=join).to_df()

    def _fetch_bins(self, range1: str | None, range2: str | None) -> pd.DataFrame | None:
        if self._is_multi_res_file() or self._f is None:
            return None

        return filter_bins(self._f.bins(), self._f.chromosomes(), range1, range2)

    def _fetch_chroms(self, range1: str | None = None, range2: str | None = None) -> Dict[str, int] | None:
        if self._f is None:
            return None

        return filter_chroms(self._f.chromosomes(), range1, range2)

    def _fetch_normalizations(self) -> List[str] | None:
        if self._is_multi_res_file():
            return hictkpy.File(self._f.path(), self._f.resolutions()[-1]).avail_normalizations()
        if self._is_single_cell_file():
            return None
        if self._f is None:
            return None

        return self._f.avail_normalizations()

    def _fetch_resolutions(self) -> List[int] | None:
        if self._is_multi_res_file():
            return self._f.resolutions()
        if self._f is None:
            return None
        return [self._f.resolution()]

    def _fetch_cells(self) -> List[str] | None:
        if self._is_single_cell_file():
            return self._f.cells()
        return None

    def _fetch_weights(self, range1: str | None, range2: str | None) -> pd.DataFrame | None:
        if not self._is_single_res_file():
            return None

        # TODO remove once hictkpy.File().weights() is available
        # https://github.com/paulsengroup/hictkpy/pull/49
        raise NotImplementedError

        data = {}
        for norm in self._f.avail_normalizations():
            data[norm] = self._f.weights(norm)

        return filter_weights(pd.DataFrame(data), self._f.chromosomes(), range1, range2)

    def dump(
        self,
        **kwargs,
    ) -> Tuple[str, pd.DataFrame | List[int] | List[str] | None]:
        table = kwargs.get("table", "pixels")
        range1 = kwargs.get("range1")
        range2 = kwargs.get("range2")

        if table == "bins":
            data = self._fetch_bins(range1, range2)
        elif table == "chroms":
            data = self._fetch_chroms(range1, range2)
        elif table == "pixels":
            balance = kwargs.get("balance", False)
            join = kwargs.get("join", False)
            cis_only = kwargs.get("cis_only", False)
            trans_only = kwargs.get("trans_only", False)
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
