#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
import time

import pandas as pd
import plotly
import plotly.express as px


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg) -> pathlib.Path:
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    cli.add_argument(
        "parquet", type=existing_file, help="Path to a .parquet file with the benchmark data to be plotted."
    )
    cli.add_argument("--force", action="store_true", default=False, help="Force overwrite existing file(s).")
    cli.add_argument(
        "--verbosity",
        type=str,
        choices={"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"},
        default="INFO",
        help="Set log verbosity.",
    )

    return cli


def generate_dummy_dataset(df: pd.DataFrame) -> pd.DataFrame:
    import numpy as np

    dfs = []
    for i in range(20):
        dff = df.copy()
        dff["tag"] = f"sha-{i:03d}"
        dff["time"] = time.time()
        dff["mean_runtime_ns"] *= np.random.normal(1.0, 0.05, 1)[0]
        dfs.append(dff)

        time.sleep(0.1)

    return pd.concat(dfs)


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())
    setup_logger(args["verbosity"])

    df = generate_dummy_dataset(pd.read_parquet(args["parquet"]))

    for key, dff in df.groupby(["test_case", "range1", "range2", "resolution"]):
        print(dff)
        print(key, dff.columns)
        fig = px.line(dff, x="tag", y="mean_runtime_ns", color="count-type")
        fig.show()
        break


if __name__ == "__main__":
    try:
        import pyarrow
    except ModuleNotFoundError as e:
        raise ImportError(
            "Please install pandas with parquet support with e.g. \"pip install 'pandas[parquet]'\""
        ) from e

    main()
