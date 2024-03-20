#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "weight1",
        type=pathlib.Path,
        help="Path to the first set of weights to be compared.",
    )
    cli.add_argument(
        "weight2",
        type=pathlib.Path,
        help="Path to the second set of weights to be compared.",
    )

    cli.add_argument(
        "weight-name",
        type=str,
        help="Name of the weights to be compared.",
    )

    cli.add_argument(
        "--tolerance",
        default=1.0e-3,
        type=float,
        help="Maximum relative error allowed.",
    )

    return cli


def compute_rel_error(x: pd.Series, y: pd.Series) -> pd.Series:
    # NaN in both Series
    mask1 = (x.isnull()) & (y.isnull())

    # NaN in only one of the Series
    mask2 = ((x.isnull()) & (~y.isnull())) | ((~x.isnull()) & (y.isnull()))

    err = ((x - y) / y).abs() * 100
    err[mask1] = 0.0
    err[mask2] = 1.0

    return err


def main() -> int:
    args = vars(make_cli().parse_args())

    name = args["weight-name"]
    tol = args["tolerance"]

    w1 = pd.read_table(args["weight1"])[name]
    w2 = pd.read_table(args["weight2"])[name]

    rerr = compute_rel_error(w1, w2)

    df = pd.DataFrame({f"{name}1": w1, f"{name}2": w2, "rel_err": rerr})
    df = df[rerr > tol]

    if len(df) != 0:
        df.to_csv(sys.stdout, index=False, sep="\t")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
