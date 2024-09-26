# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import subprocess as sp


def version(hictk: pathlib.Path) -> str:
    ver = sp.check_output([str(hictk), "--version"], encoding="utf-8").strip()
    if ver.startswith("hictk"):
        return ver.removeprefix("hictk-")

    raise RuntimeError(f'"{hictk}" does not seem to be an hictk binary')
