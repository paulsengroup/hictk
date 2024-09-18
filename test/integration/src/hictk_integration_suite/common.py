# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Tuple


def parse_uri(s: pathlib.Path | str) -> Tuple[pathlib.Path, str | None]:
    path, _, grp = str(s).partition("::")
    if grp:
        grp = str(pathlib.Path(grp).as_posix())
    else:
        grp = None

    return pathlib.Path(path), grp
