#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import yaml


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("-f", "--file", type=str)
    cli.add_argument("-q", "--quiet", action="store_true", default=False)

    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()
    try:
        if args.file:
            with open(args.file) as f:
                yaml.load(f, yaml.Loader)
        else:
            yaml.load(sys.stdin, yaml.Loader)
    except ValueError:
        print("FAILURE! String is not valid JSON", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print("OK!", file=sys.stderr)
