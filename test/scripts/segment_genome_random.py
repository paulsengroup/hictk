#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import random
import argparse
import pathlib
import pandas as pd


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "chrom-sizes",
        type=pathlib.Path,
        help="Path to a .chrom.sizes file.",
    )

    cli.add_argument(
        "bin-size",
        type=positive_int,
        help="Average bin size.",
    )

    cli.add_argument(
        "--seed",
        type=positive_int,
        default=2074288341,
        help="Seed used to initialize the PRNG.",
    )

    return cli


def read_chrom_sizes(path_to_chrom_sizes: pathlib.Path):
    return pd.read_table(path_to_chrom_sizes, names=["chrom", "length"])


def segment_chromosome(chrom: str, size: int, bin_size: int):
    start = 0
    end = min(start + int(random.gauss(bin_size, 0.1 * bin_size)), size)

    print(f"{chrom}\t{start}\t{end}")

    while end < size:
        start = end
        end = min(start + int(random.gauss(bin_size, 0.1 * bin_size)), size)
        print(f"{chrom}\t{start}\t{end}")


def main():
    args = vars(make_cli().parse_args())
    random.seed(args["seed"])

    chrom_sizes = read_chrom_sizes(args["chrom-sizes"])
    for _, (chrom, size) in chrom_sizes.iterrows():
        segment_chromosome(chrom, size, args["bin-size"])


if __name__ == "__main__":
    main()
