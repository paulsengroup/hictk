#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
import subprocess as sp
import textwrap
from typing import List


def existing_dir(arg):
    if (path := pathlib.Path(arg)).is_dir():
        return path

    raise FileNotFoundError(f'Path "{arg}" is not an existing directory')


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("--root", type=pathlib.Path, help="Path to the repository root.")

    return cli


def infer_root_dir() -> pathlib.Path:
    path = pathlib.Path(sp.check_output(["git", "rev-parse", "--show-toplevel"], encoding="utf-8").strip())

    if not path.is_dir():
        raise RuntimeError("Unable to infer repository root!")

    return path


def extract_headers(cmakelists: pathlib.Path) -> List[str]:
    include_dir = cmakelists.parent / "include"
    if not include_dir.is_dir():
        return []
    prefix = include_dir.as_posix()

    headers = []
    for f in include_dir.rglob("*.hpp"):
        f = f.as_posix().removeprefix(prefix).lstrip("/")
        if "/impl/" not in f:
            headers.append(f)

    return list(sorted(headers))


def setup_logger(level: str = "INFO"):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    root_dir = args["root"]
    if root_dir is None:
        root_dir = infer_root_dir()

    logging.info("found repository root: %s", root_dir)

    headers = []
    for cmakelists in sorted(root_dir.rglob("src/libhictk/**/CMakeLists.txt")):
        include_dir = cmakelists.parent / "include"
        if not include_dir.is_dir():
            continue

        headers.extend(extract_headers(cmakelists))

    spdx_header = textwrap.dedent(
        """
        // Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
        //
        // SPDX-License-Identifier: MIT
        """
    ).strip()

    print(spdx_header, end="\n\n")

    for h in headers:
        print(f'#include "{h}"')

    main_src = "\nint main() { std::cout << hictk::config::version::str_long() << '\\n'; }\n"

    print(main_src)


if __name__ == "__main__":
    setup_logger()
    main()
