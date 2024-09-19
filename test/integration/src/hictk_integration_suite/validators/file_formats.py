# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import pathlib
from typing import Any, Dict

import cooler

from hictk_integration_suite.common import parse_uri


def is_cooler(uri) -> bool:
    try:
        return cooler.fileops.is_cooler(pathlib.Path(uri).as_posix())
    except:  # noqa
        return False


def is_mcool(uri) -> bool:
    try:
        return cooler.fileops.is_multires_file(uri)
    except:  # noqa
        return False


def is_scool(uri) -> bool:
    try:
        return cooler.fileops.is_scool_file(uri)
    except:  # noqa
        return False


def is_hic(uri) -> bool:
    try:
        with open(uri, "rb") as f:
            magic_string = f.read(4).decode("utf-8")
        return magic_string == "HIC\0"
    except:  # noqa
        return False


def is_multires(uri) -> bool:
    return is_hic(uri) or is_mcool(uri)


def is_variable_bin_size(uri) -> bool:
    uri = pathlib.Path(uri).as_posix()
    if not is_cooler(uri):
        return False

    bin_type = cooler.Cooler(uri).info.get("bin-type", "fixed")
    return bin_type == "variable"


def get_format(uri) -> str:
    if is_cooler(uri):
        return "cool"
    if is_mcool(uri):
        return "mcool"
    if is_scool(uri):
        return "scool"
    if is_hic(uri):
        return "hic"

    return "unknown"
