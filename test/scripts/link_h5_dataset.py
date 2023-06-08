#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import pathlib

import cooler


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "cooler",
        type=pathlib.Path,
        help="Path to a Cooler file (URI syntax supported).",
    )
    cli.add_argument("src", type=str, help="Source path.")
    cli.add_argument("dest", type=str, help="Destination path.")

    return cli


def main():
    args = vars(make_cli().parse_args())

    with cooler.Cooler(str(args["cooler"])).open("r+") as f:
        f[args["dest"]] = f[args["src"]]


if __name__ == "__main__":
    main()
