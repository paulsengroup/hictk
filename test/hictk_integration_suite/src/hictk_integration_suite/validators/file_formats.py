# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import cooler


def is_cooler(uri) -> bool:
    try:
        return cooler.fileops.is_cooler(uri)
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
            magic_string = f.read(3).decode("utf-8")
        return magic_string == "HIC"
    except:  # noqa
        return False


def is_multires(uri) -> bool:
    return is_hic(uri) or is_mcool(uri)


def check_format(uri) -> str:
    if is_cooler(uri):
        return "cool"
    if is_mcool(uri):
        return "mcool"
    if is_scool(uri):
        return "scool"
    if is_hic(uri):
        return "hic"

    return "unknown"
