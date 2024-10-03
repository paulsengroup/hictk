# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
import subprocess as sp


def version(hictk: pathlib.Path) -> str:
    env_variables = os.environ.copy()
    if "LLVM_PROFILE_FILE" in env_variables:
        env_variables["LLVM_PROFILE_FILE"] = env_variables["LLVM_PROFILE_FILE"].replace("%id", "%p")

    ver = sp.check_output([str(hictk), "--version"], encoding="utf-8", env=env_variables).strip()
    if ver.startswith("hictk"):
        return ver.removeprefix("hictk-")

    raise RuntimeError(f'"{hictk}" does not seem to be an hictk binary')
