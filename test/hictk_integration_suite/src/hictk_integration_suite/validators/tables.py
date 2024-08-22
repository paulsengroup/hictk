# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from typing import Dict, List, Sized

import numpy as np
import pandas as pd


def _compare_lists(expected: List, found: List, name: str, sort: bool = True) -> Dict[str, str]:
    if len(expected) != 0 and len(found) == 0:
        return {f"no {name} were fetched": ""}

    if sort:
        expected.sort()
        found.sort()

    if expected != found:
        return {f"found differences in {name}": f"expected {expected}, found {found}"}

    return {}


def _empty_sized_to_none(data: Sized | None) -> Sized | None:
    try:
        len(data)
        return data
    except:  # noqa
        return None


def _handle_common_errors(expected: Sized | None, found: Sized | None, title: str) -> Dict[str, str]:
    expected = _empty_sized_to_none(expected)
    found = _empty_sized_to_none(found)

    if expected is None and found is None:
        return {}

    if expected is not None and found is None:
        return {f"no {title} were fetched": ""}

    if expected is None and found is not None:
        return {"unexpected output to stdout": f"expected no output, read {len(found)} {title}"}

    return {}


def compare_bins(
    expected: pd.DataFrame | None,
    found: pd.DataFrame | None,
    ignore_index: bool = True,
    sort_before_compare: bool = False,
) -> Dict[str, str]:
    errors = _handle_common_errors(expected, found, "bins")
    if len(errors) != 0:
        return errors

    if expected.columns.tolist() != found.columns.tolist():
        return {"column name mismatch": f"expected {expected.columns.tolist()}, found {found.columns.tolist()}"}

    columns = ["chrom", "start", "end"]
    if sort_before_compare:
        expected = expected.sort_values(columns)
        found = found.sort_values(columns)

    coord_differences = pd.concat([expected, found]).drop_duplicates(
        subset=columns, keep=False, ignore_index=ignore_index
    )

    if len(coord_differences) != 0:
        return {"found differences in pixel coordinates": f"found {len(coord_differences)} differences"}

    return {}


def compare_chroms(expected: Dict[str, int] | None, found: pd.DataFrame | None) -> Dict[str, str]:
    errors = _handle_common_errors(expected, found, "chromosomes")
    if len(errors) != 0:
        return errors

    found = found.set_index("chrom")["size"].to_dict()
    if list(expected.keys()) != list(found.keys()):
        return {"chromosome names mismatch": f"expected {list(expected.keys())}, found {list(found.keys())}"}

    if expected != found:
        return {"chromosome sizes mismatch": f"expected {expected}, found {found}"}

    return {}


def compare_pixels(
    expected: pd.DataFrame,
    found: pd.DataFrame,
    rtol: float = 1.0e-5,
    ignore_index: bool = True,
    sort_before_compare: bool = False,
) -> Dict[str, str]:
    assert 0 <= rtol <= 1.0
    errors = _handle_common_errors(expected, found, "pixels")
    if len(errors) != 0:
        return errors

    if expected.columns.tolist() != found.columns.tolist():
        return {"column name mismatch": f"expected {expected.columns.tolist()}, found {found.columns.tolist()}"}

    if (expected.dtypes != found.dtypes).any():
        return {"column data type mismatch": f"expected {expected.dtypes}, found {found.dtypes}"}

    if len(expected) != len(found):
        return {"record number mismatch": f"expected {len(expected)} records, found {len(found)}"}

    columns = expected.columns.tolist()
    columns.remove("count")
    if sort_before_compare:
        expected = expected.sort_values(columns)
        found = found.sort_values(columns)

    expected = expected.copy()
    found = found.copy()
    expected["type"] = "expected"
    found["type"] = "found"

    coord_differences = pd.concat([expected, found]).drop_duplicates(
        subset=columns, keep=False, ignore_index=ignore_index
    )

    if len(coord_differences) != 0:
        return {"found differences in pixel coordinates": f"found {len(coord_differences)} differences"}

    count_differences = (~np.isclose(expected["count"], found["count"], rtol=rtol, equal_nan=True)).sum()
    if count_differences != 0:
        return {"found differences in pixel counts": f"found {count_differences} differences"}

    return {}


def compare_normalizations(expected: List[str], found: pd.DataFrame) -> Dict[str, str]:
    errors = _handle_common_errors(expected, found, "normalizations")
    if len(errors) != 0:
        return errors
    return _compare_lists(expected, found["normalization"].tolist(), "normalizations")


def compare_resolutions(expected: List[int], found: pd.DataFrame) -> Dict[str, str]:
    errors = _handle_common_errors(expected, found, "resolutions")
    if len(errors) != 0:
        return errors
    return _compare_lists(expected, found["resolution"].tolist(), "resolutions")


def compare_cells(expected: List[str], found: pd.DataFrame) -> Dict[str, str]:
    errors = _handle_common_errors(expected, found, "cells")
    if len(errors) != 0:
        return errors
    return _compare_lists(expected, found["cell"].tolist(), "cells")


def compare_weights(
    expected: pd.DataFrame,
    found: pd.DataFrame,
    rtol: float = 1.0e-5,
    sort_before_compare: bool = True,
) -> Dict[str, str]:
    assert 0 <= rtol <= 1.0
    errors = _handle_common_errors(expected, found, "weights")
    if len(errors) != 0:
        return errors

    if sort_before_compare:
        columns = list(sorted(expected.columns.tolist()))
        expected = expected[columns]
        columns = list(sorted(found.columns.tolist()))
        found = found[columns]

    if expected.columns.tolist() != found.columns.tolist():
        return {"column name mismatch": f"expected {expected.columns.tolist()}, found {found.columns.tolist()}"}

    expected = expected.astype(float)
    found = found.astype(float)

    if len(expected) != len(found):
        return {"record number mismatch": f"expected {len(expected)} records, found {len(found)}"}

    errors = {}
    for col in expected:
        differences = (~np.isclose(expected[col], found[col], rtol=rtol, equal_nan=True)).sum()
        if differences != 0:
            errors |= {f'found differences in "{col}" weights': f"found {differences} differences"}

    return errors
