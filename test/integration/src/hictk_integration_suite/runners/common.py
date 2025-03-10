# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
from typing import Dict

import cooler
import numpy as np
import pandas as pd


def normalize_df_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "chrom1" in df.columns:
        df[["chrom1", "chrom2"]] = df[["chrom1", "chrom2"]].astype(str)
    if "chrom" in df.columns:
        df["chrom"] = df["chrom"].astype(str)

    int_types = (int, np.uint8, np.uint16, np.uint32, np.uint64, np.int8, np.int16, np.int32, np.int64)
    float_types = (float, np.float16, np.float32, np.float64, np.longdouble)
    columns = df.select_dtypes(include=int_types).columns.tolist()
    df[columns] = df[columns].astype(np.int64)

    columns = df.select_dtypes(include=float_types).columns.tolist()
    df[columns] = df[columns].astype(np.float64)

    return df


def filter_bins(df: pd.DataFrame, chroms: Dict[str, int], range1: str | None, range2: str | None) -> pd.DataFrame:
    df = normalize_df_dtypes(df)
    if range1 is None:
        assert range2 is None
        return df

    chrom1, start1, end1 = cooler.api.parse_region(range1, chroms)
    df1 = df[(df["chrom"] == chrom1) & (df["start"] >= start1) & (df["start"] < end1)]
    if range2 is None or range1 == range2:
        return df1

    chrom2, start2, end2 = cooler.api.parse_region(range2, chroms)
    df2 = df[(df["chrom"] == chrom2) & (df["start"] >= start2) & (df["start"] < end2)]

    return pd.concat([df1, df2])


def filter_chroms(chroms: Dict[str, int], range1: str | None, range2: str | None) -> Dict[str, int]:
    if range1 is None:
        assert range2 is None
        return chroms

    chrom1, _, _ = range1.partition(":")
    if range2 is None:
        return {chrom1: chroms[chrom1]}

    chrom2, _, _ = range2.partition(":")
    return {chrom1: chroms[chrom1], chrom2: chroms[chrom2]}


def filter_weights(df: pd.DataFrame, chroms: Dict[str, int], range1: str | None, range2: str | None) -> pd.DataFrame:
    df = normalize_df_dtypes(df)

    if range1 is None:
        assert range2 is None
        return df

    chrom1, start1, end1 = cooler.api.parse_region(range1, chroms)
    df1 = df[(df["chrom"] == chrom1) & (df["start"] >= start1) & (df["start"] < end1)]
    if range2 is None or range1 == range2:
        return df1

    chrom2, start2, end2 = cooler.api.parse_region(range2, chroms)
    df2 = df[(df["chrom"] == chrom2) & (df["start"] >= start2) & (df["start"] < end2)]

    return pd.concat([df1, df2])
