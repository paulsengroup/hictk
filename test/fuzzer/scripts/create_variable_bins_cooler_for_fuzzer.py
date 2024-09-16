#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import gzip
import logging
import multiprocessing as mp
import pathlib
import random
import re
import subprocess as sp
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor as Pool
from typing import Dict

import cooler
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    cli.add_argument(
        "pairs",
        type=existing_file,
        help="Path to a .pairs file to use as input.",
    )
    cli.add_argument(
        "output-prefix",
        type=pathlib.Path,
        help="Path prefix (including parent folder(s)) to use for output.",
    )
    cli.add_argument(
        "--chrom-sizes",
        type=existing_file,
        help="Path to a .chrom.sizes file with the list of chromosomes to use as reference.\n"
        "When not provided, chromosome sizes are extracted from the .pairs file header.",
    )
    cli.add_argument(
        "--bin-size-avg",
        type=float,
        default=50_000,
        help="Mean value used to sample bin intervals from a normal distribution.",
    )
    cli.add_argument(
        "--bin-size-std",
        type=float,
        default=50_000,
        help="Standard deviation value used to sample bin intervals from a normal distribution.",
    )
    cli.add_argument(
        "--tmpdir",
        type=pathlib.Path,
        default=pathlib.Path(tempfile.gettempdir()),
        help="Path to a folder to use for temporary file(s).",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    cli.add_argument(
        "--seed",
        type=int,
        default=1222139576,
        help="Seed used to initialize the PRNG.",
    )
    cli.add_argument(
        "--verbosity",
        type=str,
        choices={"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"},
        default="INFO",
        help="Set log verbosity.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def ingest_pairs(
    pairs: pathlib.Path,
    bin_table: pathlib.Path,
    out_prefix: pathlib.Path,
    tmpdir: pathlib.Path,
    force: bool,
) -> pathlib.Path:
    dest = pathlib.Path(f"{out_prefix}.variable.cool")

    if force:
        dest.unlink(missing_ok=True)

    if dest.exists():
        raise RuntimeError(
            f'refusing to overwrite existing file "{dest}". Pass --force to overwrite existing file(s).'
        )

    cmd = [
        "cooler",
        "cload",
        "pairs",
        bin_table,
        pairs,
        dest,
        "-c1",
        "2",
        "-p1",
        "3",
        "-c2",
        "4",
        "-p2",
        "5",
        "--temp-dir",
        tmpdir,
    ]
    sp.check_call(cmd),
    return dest


def file_is_gzipped(path: pathlib.Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except:  # noqa
        return False


def try_open_gzipped_file(path: pathlib.Path):
    if file_is_gzipped(path):
        return gzip.open(path, "rt")

    return open(path, "r")


def read_chrom_sizes_from_header(path_to_pairs: pathlib.Path) -> Dict[str, int]:
    chroms = {}
    pattern = re.compile(r"#chromsize:\s+([\w_\-]+)\s+(\d+)$")
    with try_open_gzipped_file(path_to_pairs) as f:
        for line in f:
            if not line.startswith("#"):
                break
            matches = pattern.search(line)
            if matches:
                chroms[matches.group(1)] = int(matches.group(2))

    if len(chroms) == 0:
        raise RuntimeError(f"Unable to read any chromosome size from {path_to_pairs}")
    return chroms


def generate_bin_size(avg: float, std: float) -> int:
    bin_size = 0
    while bin_size < 1:
        bin_size = random.gauss(avg, std)

    return int(bin_size)


def generate_bin_table(
    chrom_sizes: Dict[str, int], bin_size_avg: float, bin_size_std: float
) -> pd.DataFrame:
    chroms = []
    starts = []
    ends = []

    for chrom, size in chrom_sizes.items():
        starts_ = []
        ends_ = []

        while True:
            if len(ends_) == 0:
                bin_start = 0
            else:
                bin_start = ends_[-1]

            bin_size = generate_bin_size(bin_size_avg, bin_size_std)
            bin_end = min(bin_start + bin_size, size)

            starts_.append(bin_start)
            ends_.append(bin_end)

            if bin_end == size:
                break

        chroms.extend([chrom] * len(starts_))
        starts.extend(starts_)
        ends.extend(ends_)

    return pd.DataFrame({"chrom": chroms, "start": starts, "end": ends})


def main():
    args = vars(make_cli().parse_args())
    setup_logger(args["verbosity"])

    out_prefix = args["output-prefix"]
    pairs = args["pairs"]
    nproc = args["nproc"]
    force = args["force"]

    with tempfile.TemporaryDirectory(
        prefix=str(args["tmpdir"] / "hictk-")
    ) as tmpdir, Pool(
        max_workers=nproc, initializer=setup_logger, initargs=(args["verbosity"],)
    ) as pool:
        tmpdir = pathlib.Path(tmpdir)
        out_prefix.parent.mkdir(parents=True, exist_ok=True)

        if args["chrom_sizes"] is None:
            chrom_sizes = read_chrom_sizes_from_header(pairs)
        else:
            chrom_sizes = (
                pd.read_table(args["chrom_sizes"], names=["chrom", "length"])
                .set_index("chrom")
                .to_dict()
            )

        random.seed(args["seed"])
        bin_table = generate_bin_table(
            chrom_sizes, args["bin_size_avg"], args["bin_size_std"]
        )

        bin_table_bed = tmpdir / "bins.bed"
        bin_table.to_csv(bin_table_bed, sep="\t", index=False, header=False)

        # pairs to cool
        clr = ingest_pairs(
            pairs,
            bin_table_bed,
            out_prefix,
            tmpdir,
            force,
        )

        logging.info('balancing cooler at "%s"...', clr)
        t0 = time.time()
        cooler.balance_cooler(cooler.Cooler(str(clr)), store=True, map=pool.map)
        t1 = time.time()
        logging.info('DONE! balancing cooler at "%s" took %ss', clr, t1 - t0)


if __name__ == "__main__":
    main()
