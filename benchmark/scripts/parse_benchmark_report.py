#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import json
import logging
import pathlib
import sys
import xml.etree.ElementTree as ET
from typing import IO, Any, Dict, List

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def existing_file(arg) -> pathlib.Path | str:
        if arg == "-":
            return "stdin"
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    cli.add_argument(
        "xml",
        type=existing_file,
        help='Path to the benchmark report produced by Catch2 in XML format.\nPass "-" to read data from stdin.',
    )
    cli.add_argument("parquet", type=pathlib.Path, help="Path where to store the output in .parquet format.")
    cli.add_argument("--force", action="store_true", default=False, help="Force overwrite existing file(s).")
    cli.add_argument(
        "--verbosity",
        type=str,
        choices={"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"},
        default="INFO",
        help="Set log verbosity.",
    )

    return cli


def parse_mean_node(data: Dict[str, str]) -> float:
    return int(float(data["value"]))


def parse_stddev_node(data: Dict[str, str]) -> float:
    return int(float(data["value"]))


def parse_benchmark_result(data: Dict[str, str]) -> Dict[str, Any]:
    name = json.loads(data["name"])
    return name | {"clock_resolution_ns": float(data["clockResolution"])}


def parse_test_case(data: Dict[str, str]) -> Dict[str, Any]:
    return json.loads(data["name"])


def log_event_start(element):
    if element.tag == "TestCase":
        logging.debug('begin parsing TestCase "%s"...', element.attrib.get("name"))
    if element.tag == "BenchmarkResults":
        try:
            name = json.loads(element.attrib["name"])
        except:  # noqa
            name = element.attrib["name"]
        logging.debug('begin parsing BenchmarkResults "%s"...', name)


def finalize_test_case(data: Dict[str, str], records: List[Dict]) -> List[Dict]:
    test_case = parse_test_case(data)
    return [test_case | r for r in records]


def finalize_benchmark_results(data: Dict[str, str], record: Dict[str, Any]) -> Dict[str, Any]:
    result = parse_benchmark_result(data)
    return result | record


def parse_document(source: pathlib.Path | IO) -> pd.DataFrame:
    records = []
    test_case_records = []
    benchmark_result = {}
    for event, element in ET.iterparse(source, ["start", "end"]):
        if event == "start":
            log_event_start(element)
        else:
            assert event == "end"
            if element.tag == "TestCase":
                records.extend(finalize_test_case(element.attrib, test_case_records))
                test_case_records = []
                benchmark_result = {}
            elif element.tag == "BenchmarkResults":
                test_case_records.append(finalize_benchmark_results(element.attrib, benchmark_result))
                benchmark_result = {}
            elif element.tag == "mean":
                benchmark_result["mean_runtime_ns"] = parse_mean_node(element.attrib)
            elif element.tag == "standardDeviation":
                benchmark_result["stddev_runtime_ns"] = parse_stddev_node(element.attrib)

    if len(benchmark_result) != 0:
        raise RuntimeError("failed to parse BenchmarkResults node: end token is missing!")

    if len(test_case_records) != 0:
        raise RuntimeError("failed to parse TestCase node: end token is missing!")

    if len(records) == 0:
        raise RuntimeError("no records were parsed. Is the document empty?")

    df = pd.DataFrame(records)

    for col in df.filter(regex=r".*_ns$"):
        df[col] = pd.to_timedelta(df[col], unit="ns")
        df = df.rename({col: col.removesuffix("_ns")})

    return df


def raise_file_exist_except(path: pathlib.Path):
    raise RuntimeError(f'refusing to overwrite file "{path}". Pass --force to overwrite existing file(s).')


def setup_logger(level):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())
    setup_logger(args["verbosity"])

    src = args["xml"]
    dest = args["parquet"]

    if dest.exists() and not args["force"]:
        raise_file_exist_except(dest)

    dest.unlink(missing_ok=True)

    if src == "stdin":
        df = parse_document(sys.stdin)
    else:
        df = parse_document(src)

    logging.info(f'successfully parsed {len(df)} records from "{src}"!')

    if dest.exists():
        raise_file_exist_except(dest)

    logging.info(f'writing records to file "{dest}"...')
    dest.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(dest)


if __name__ == "__main__":
    try:
        import pyarrow
    except ModuleNotFoundError as e:
        raise ImportError(
            "Please install pandas with parquet support with e.g. \"pip install 'pandas[parquet]'\""
        ) from e

    main()
