#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import re


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser("Patch hictkpy's CMakeLists.txt.")

    cli.add_argument(
        "--cmakelists",
        type=str,
        required=True,
        help="Path to the CMakeLists.txt to be patched.",
    )
    cli.add_argument(
        "--git-tag",
        type=str,
        required=True,
        help="Git tag to be used to override the version of hictk used by FetchContent.",
    )
    cli.add_argument(
        "-i",
        "--inplace",
        action="store_true",
        default=False,
        help="Apply the patch to the input CMakeLists.txt.",
    )

    return cli


if __name__ == "__main__":
    args = vars(make_cli().parse_args())
    repo_url = "https://github.com/paulsengroup/hictk.git"
    git_tag = args["git_tag"]

    pattern = re.compile(
        r"(FetchContent_Declare\(\s+hictk\s+)URL.*(\s+)URL_HASH.*(\s+EXCLUDE_FROM_ALL\s+SYSTEM\))"
    )

    with open(args["cmakelists"]) as f:
        data = "".join(f.readlines())

    data = pattern.sub(f"\\1GIT_REPOSITORY {repo_url}\\2GIT_TAG {git_tag}\\3", data)

    if args["inplace"]:
        with open(args["cmakelists"], "w") as f:
            f.write(data)
    else:
        print(data)
