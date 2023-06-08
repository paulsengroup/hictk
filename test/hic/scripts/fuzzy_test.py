#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import ctypes
import gc
import itertools
import logging
import multiprocessing as mp
import pathlib
import random
import shlex
import shutil
import subprocess as sp
import sys
import time
from typing import Dict, Tuple, Union

import hicstraw
import numpy as np
import pandas as pd


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    def valid_fraction(arg):
        if (n := float(arg)) >= 0 and n <= 1:
            return n

        raise ValueError("Not a number between 0 and 1")

    def valid_executable(arg):
        if (cmd := shutil.which(arg)) is not None:
            return pathlib.Path(cmd)

        raise FileNotFoundError(f'Unable to find executable "{arg}"')

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "hic",
        type=pathlib.Path,
        help="Path to a .hic file.",
    )

    cli.add_argument("--resolution", type=positive_int, help="Matrix resolution.")

    cli.add_argument(
        "--1d-to-2d-query-ratio",
        type=valid_fraction,
        default=0.33,
        help="Ratio of 1D to 2D queries. Use 0 or 1 to only test 1D or 2D queries.",
    )

    cli.add_argument("--duration", type=positive_int, default=60, help="Duration in seconds.")
    cli.add_argument(
        "--path-to-hicxx-dump",
        type=valid_executable,
        default="hicxx_dump",
        help="Path to hicxx_example binary.",
    )
    cli.add_argument("--query-length-avg", type=float, default=5_000_000, help="Average query size.")
    cli.add_argument(
        "--query-length-std",
        type=float,
        default=1_000_000,
        help="Standard deviation for query size.",
    )
    cli.add_argument("--balance", type=str, default="NONE", help="Name of the dataset to use for balancing.")
    cli.add_argument("--seed", type=int, default=2074288341)
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    return cli


def parse_ucsc_string(s: str) -> Tuple[str, int, int]:
    chrom, start, end = s.replace("-", ":").split(":")
    return chrom, int(start), int(end)


def hic_records_to_df(chrom1: str, chrom2: str, resolution: int, chromosomes, records) -> pd.DataFrame:
    pos1 = []
    pos2 = []
    counts = []

    size1 = chromosomes[chrom1]
    size2 = chromosomes[chrom2]

    for r in records:
        pos1.append(r.binX)
        pos2.append(r.binY)
        counts.append(r.counts)

    df = pd.DataFrame(
        {
            "chrom1": [chrom1] * len(counts),
            "start1": pos1,
            "end1": 0,
            "chrom2": [chrom2] * len(counts),
            "start2": pos2,
            "end2": 0,
            "count": counts,
        }
    )

    df["end1"] = np.minimum(df["start1"] + resolution, size1)
    df["end2"] = np.minimum(df["start2"] + resolution, size2)
    return df


def hicstraw_dump(
    selectors, chromosomes, resolution, chrom1: str, start1: int, end1: int, chrom2: str, start2: int, end2: int
) -> pd.DataFrame:
    logging.debug("[hicstraw] running query for %s:%d-%d, %s:%d-%d...", chrom1, start1, end1, chrom2, start2, end2)
    sel = selectors[(chrom1, chrom2)]

    # Fix up queries. hicstraw doesn't include partially overlapping bins on the left side
    start11 = (start1 // resolution) * resolution
    start22 = (start2 // resolution) * resolution

    df = hic_records_to_df(chrom1, chrom2, resolution, chromosomes, sel.getRecords(start11, end1, start22, end2))

    return df[(df["start1"] < end1) & (df["start2"] < end2)].set_index(
        ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    )


def hicxx_dump(
    hicxx_dump_bin: pathlib.Path,
    path_to_hic_file: pathlib.Path,
    balance: Union[bool, str],
    resolution: int,
    chrom1: str,
    start1: int,
    end1: int,
    chrom2: str,
    start2: int,
    end2: int,
) -> pd.DataFrame:
    cmd = [
        shutil.which(str(hicxx_dump_bin)),
        "observed",
        balance,
        str(path_to_hic_file),
        f"{chrom1}:{start1}-{end1}",
        f"{chrom2}:{start2}-{end2}",
        "BP",
        str(resolution),
    ]

    cmd = shlex.split(" ".join(str(tok) for tok in cmd))
    logging.debug("[hicxx] Running %s...", cmd)
    sp.check_output(cmd, stderr=sp.STDOUT)

    with sp.Popen(cmd, stdin=None, stderr=sp.PIPE, stdout=sp.PIPE) as hicxx_dump:
        df = pd.read_table(
            hicxx_dump.stdout,
            names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"],
        )
        hicxx_dump.communicate()
        if (code := hicxx_dump.returncode) != 0:
            print(hicxx_dump.stderr, file=sys.stderr)
            raise RuntimeError(f"{cmd} terminated with code {code}")

    df["count"] = df["count"].astype(float)
    return df.set_index(["chrom1", "start1", "end1", "chrom2", "start2", "end2"])


def read_chrom_sizes(path_to_hic_file: pathlib.Path) -> Dict[str, int]:
    chroms = {}
    for chrom in hicstraw.HiCFile(str(path_to_hic_file)).getChromosomes():
        if str(chrom.name).lower() != "all":
            chroms[chrom.name] = chrom.length
    return chroms


def generate_query_1d(chroms, weights: np.ndarray, mean_length: float, stddev_length: float) -> Tuple[str, int, int]:
    chrom_name, chrom_size = random.choices(chroms, weights=weights, k=1)[0]

    query_length = max(2.0, random.gauss(mu=mean_length, sigma=stddev_length))

    center_pos = random.randint(0, chrom_size)
    start_pos = max(0.0, center_pos - (query_length / 2))
    end_pos = min(chrom_size, start_pos + query_length)

    return chrom_name, int(start_pos), int(end_pos)


def generate_query_2d(
    chroms,
    weights: np.ndarray,
    ranks: Dict[str, int],
    mean_length: float,
    stddev_length: float,
) -> Tuple[Tuple[str, int, int], Tuple[str, int, int]]:
    q1 = generate_query_1d(chroms, weights, mean_length, stddev_length)
    q2 = generate_query_1d(chroms, weights, mean_length, stddev_length)

    chrom1, start1, coord1 = q1
    chrom2, start2, coord2 = q2

    if ranks[chrom1] > ranks[chrom2]:
        q1, q2 = q2, q1

    if chrom1 == chrom2:
        if int(start1) > int(start2):
            q1, q2 = q2, q1

    return q1, q2


def find_differences(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    df = df1.merge(df2, how="outer", left_index=True, right_index=True, suffixes=("1", "2"), sort=False)
    # We're mapping False to None so that we can more easily drop identical rows with dropna()
    df["count_close_enough"] = pd.Series(np.isclose(df["count1"], df["count2"])).map({False: None})

    # We're dropping the counts to avoid incorrectly flagging rows with nan as counts
    return df.drop(columns=["count1", "count2"]).dropna()


def results_are_identical(worker_id, q1, q2, expected, found) -> bool:
    if len(expected) != len(found):
        logging.warning(
            "[%d] %s, %s: FAIL! Expected %d nnz, found %d!",
            worker_id,
            q1,
            q2,
            len(expected),
            len(found),
        )
        return False

    if len(expected) != 0:
        diff = find_differences(expected, found)
        if len(diff) != 0:
            logging.warning(
                "[%d] %s, %s (%d nnz): FAIL! Found %d differences!",
                worker_id,
                q1,
                q2,
                len(expected),
                len(diff),
            )
            return False

    logging.debug("[%d] %s, %s (%d nnz): OK!", worker_id, q1, q2, len(expected))
    return True


def seed_prng(worker_id: int, seed):
    seed = hash(tuple([worker_id, seed]))
    logging.info("[%d] seed: %d", worker_id, seed)
    random.seed(seed)


def init_matrix_selectors(
    path_to_hic: pathlib.Path, norm: str, resolution: int
) -> Dict[Tuple[str, str], hicstraw.MatrixZoomData]:
    f = hicstraw.HiCFile(str(path_to_hic))

    selectors = {}
    for chrom1, chrom2 in itertools.product(f.getChromosomes(), repeat=2):
        if str(chrom1.name).lower() == "all" or str(chrom2.name).lower() == "all":
            continue
        if chrom1.index <= chrom2.index:
            sel = f.getMatrixZoomData(str(chrom1.name), str(chrom2.name), "observed", norm, "BP", resolution)
            selectors[(chrom1.name, chrom2.name)] = sel

    return selectors


def worker(
    path_to_hic: pathlib.Path,
    path_to_hicxx_dump: pathlib.Path,
    resolution: int,
    chroms_flat,
    chrom_ranks,
    query_length_mu: float,
    query_length_std: float,
    _1d_to_2d_query_ratio: float,
    balance: str,
    seed: int,
    worker_id: int,
    end_time,
) -> Tuple[int, int]:
    global early_return

    num_failures = 0
    num_queries = 0

    if balance is None or balance == "raw":
        balance = False

    try:
        seed_prng(worker_id, seed)

        chrom_sizes = np.array([n for _, n in chroms_flat], dtype=int)
        weights = chrom_sizes / chrom_sizes.sum()

        selectors = None
        chromosomes = read_chrom_sizes(path_to_hic)

        for i in itertools.count():
            if time.time() >= end_time:
                break

            if i % 50000 == 0:
                del selectors
                gc.collect()
                selectors = init_matrix_selectors(path_to_hic, balance, resolution)

            if early_return.value:
                logging.debug(
                    "[%d] early return signal received. Returning immediately!",
                    worker_id,
                )
                break

            if _1d_to_2d_query_ratio <= random.random():
                q1, q2 = generate_query_2d(
                    chroms_flat,
                    weights,
                    chrom_ranks,
                    mean_length=query_length_mu,
                    stddev_length=query_length_std,
                )
            else:
                q1 = generate_query_1d(
                    chroms_flat,
                    weights,
                    mean_length=query_length_mu,
                    stddev_length=query_length_std,
                )
                q2 = q1

            num_queries += 1
            expected = hicstraw_dump(selectors, chromosomes, resolution, *q1, *q2)
            found = hicxx_dump(path_to_hicxx_dump, path_to_hic, balance, resolution, *q1, *q2)

            if not results_are_identical(worker_id, q1, q2, expected, found):
                num_failures += 1

    except:
        logging.debug(
            "[%d] exception raised in worker process. Sending early return signal!",
            worker_id,
        )
        early_return.value = True
        raise

    return num_queries, num_failures


def main():
    args = vars(make_cli().parse_args())

    chroms = read_chrom_sizes(args["hic"])
    chrom_ranks = {chrom: i for i, chrom in enumerate(chroms.keys())}
    chroms_flat = list(chroms.items())

    end_time = time.time() + args["duration"]

    with mp.Pool(args["nproc"]) as pool:
        results = pool.starmap(
            worker,
            zip(
                itertools.repeat(args["hic"]),
                itertools.repeat(args["path_to_hicxx_dump"]),
                itertools.repeat(args["resolution"]),
                itertools.repeat(chroms_flat),
                itertools.repeat(chrom_ranks),
                itertools.repeat(args["query_length_avg"]),
                itertools.repeat(args["query_length_std"]),
                itertools.repeat(args["1d_to_2d_query_ratio"]),
                itertools.repeat(args["balance"]),
                itertools.repeat(args["seed"]),
                range(args["nproc"]),
                itertools.repeat(end_time),
            ),
            chunksize=1,
        )

    num_queries = sum((n for n, _ in results))
    num_failures = sum((n for _, n in results))
    num_passes = num_queries - num_failures
    if num_failures == 0:
        lvl = logging.INFO
    else:
        lvl = logging.WARN

    logging.log(
        lvl,
        "Score: %.4g%% (%d successes and %d failures).",
        100 * num_passes / num_queries,
        num_passes,
        num_failures,
    )

    return num_failures != 0


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    early_return = mp.Value(ctypes.c_bool, False)
    sys.exit(main())
