#!/usr/bin/env python3

# Copyright (c) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import re
import shutil
import subprocess as sp
import textwrap
from typing import Tuple


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def valid_executable(s: str) -> pathlib.Path:
        if shutil.which(s):
            return pathlib.Path(s)

        if s == "hictk":
            raise argparse.ArgumentTypeError("Unable to find hictk in your PATH.")

        raise argparse.ArgumentTypeError(f'"{s}" is not a valid executable.')

    cli.add_argument(
        "--hictk",
        type=valid_executable,
        default=pathlib.Path("hictk"),
        required=False,
        help="Path to hictk's executable.",
    )

    return cli


def generate_main_header():
    header = """
    ..
       Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
       SPDX-License-Identifier: MIT

    CLI Reference
    #############

    For an up-to-date list of subcommands and CLI options refer to ``hictk --help``.

    Subcommands
    -----------

    .. code-block:: text

    """

    print(textwrap.dedent(header))


def generate_subcommand_header(subcommand: Tuple[str]):
    subcommand = "hictk " + " ".join(subcommand)
    separator = "-" * len(subcommand)
    header = f"""

    {subcommand}
    {separator}

    .. code-block:: text
    """

    print(textwrap.dedent(header))


def sanitize(msg: str, hictk: pathlib.Path) -> str:
    msg = msg.replace(str(hictk), "hictk")
    msg = re.sub(r"^\s+$", "", msg, flags=re.MULTILINE)
    return re.sub(r"\s+$", "", msg, flags=re.MULTILINE)


def generate_main_help_msg(hictk: pathlib.Path):
    msg = sp.check_output([hictk, "--help"]).decode("utf-8")
    msg = sanitize(msg, hictk)
    print(textwrap.indent(msg, "  "))


def generate_subcommand_help_msg(hictk: pathlib.Path, subcommand: Tuple[str]):
    msg = sp.check_output([hictk, *subcommand, "--help"]).decode("utf-8")
    msg = sanitize(msg, hictk)
    print(textwrap.indent(msg, "  "))


def main():
    args = vars(make_cli().parse_args())

    hictk = args["hictk"]
    if not shutil.which(hictk):
        if hictk == "hictk":
            raise RuntimeError("Unable to find hictk in your PATH.")

        raise RuntimeError(f'"{hictk}" is not a valid executable.')

    subcommands = (
        ("balance",),
        ("balance", "ice"),
        ("balance", "scale"),
        ("balance", "vc"),
        ("convert",),
        ("dump",),
        ("fix-mcool",),
        ("load",),
        ("merge",),
        ("metadata",),
        ("rename-chromosomes",),
        ("validate",),
        ("zoomify",),
    )

    generate_main_header()
    generate_main_help_msg(hictk)

    for subcommand in subcommands:
        generate_subcommand_header(subcommand)
        generate_subcommand_help_msg(hictk, subcommand)


if __name__ == "__main__":
    main()
