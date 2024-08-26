#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import datetime
import importlib
import json
import logging
import os
import pathlib
import platform
import sys
import time
import tomllib
from typing import Any, Dict, List

import click

from hictk_integration_suite.cli.common import WorkingDirectory
from hictk_integration_suite.runners.hictk.common import version


def nproc() -> int:
    return len(os.sched_getaffinity(0))


def get_test_names(include_all: bool = True) -> List[str]:
    if include_all:
        names = ["all"]
    else:
        names = []
    names += ["balance", "convert", "dump", "metadata", "zoomify"]
    return names


def update_uris(config: Dict, data_dir: pathlib.Path) -> Dict:
    def _update_uri(uri: pathlib.Path) -> str:
        uri = str(data_dir / uri)
        path, _, grp = uri.partition("::")
        path = pathlib.Path(path)
        if not path.exists():
            raise RuntimeError(f'file "{path}" does not exists')

        path = path.resolve()
        if not grp:
            return str(path)
        return f"{path}::{grp}"

    if "files" not in config and "test-cases" not in config:
        return config

    new_config = config.copy()
    for i, mappings in enumerate(config.get("files", [])):
        for key in mappings:
            if key.endswith("uri") or key.endswith("path"):
                new_config["files"][i][key] = _update_uri(mappings[key])

    for i, mappings in enumerate(config.get("test-cases", [])):
        for key in mappings:
            if key.endswith("uri") or key.endswith("path"):
                new_config["test-cases"][i][key] = _update_uri(mappings[key])

    return new_config


def stage_input_files(config: Dict[str, Any], wd: WorkingDirectory):
    for mappings in config.get("files", []):
        for key, value in mappings.items():
            if key.endswith("uri") or key.endswith("path"):
                wd.stage_file(value, exists_ok=True)


def import_config_and_stage_files(
    path: pathlib.Path,
    data_dir: pathlib.Path | None,
    wd: WorkingDirectory,
    command: str,
) -> Dict[str, Any]:
    with open(path, "rb") as f:
        config = tomllib.load(f)

    config = config[command]
    if data_dir:
        config = update_uris(config, data_dir)

    stage_input_files(config, wd)
    return config


def parse_log_lvl(lvl: str):
    if lvl == "debug":
        return logging.DEBUG
    if lvl == "info":
        return logging.INFO
    if lvl == "warning":
        return logging.WARNING
    if lvl == "error":
        return logging.ERROR
    if lvl == "critical":
        return logging.CRITICAL
    raise NotImplementedError


def init_results(hictk_bin: pathlib.Path) -> Dict:
    res = {
        "platform": platform.platform(),
        "arch": platform.machine(),
        "hictk-version": version(hictk_bin),
        "date": datetime.datetime.now().isoformat(),
        "results": {
            "pass": 0,
            "fail": 0,
        },
    }

    return res


def parse_test_suites(s: str) -> List[str]:
    suites = s.split(",")
    if isinstance(suites, str):
        suites = [suites]
    else:
        suites = set(suites)

    if "all" in suites:
        suites = get_test_names(include_all=False)
    else:
        avail_suites = set(get_test_names(include_all=True))
        for s in suites:
            if s not in avail_suites:
                valid_suites = ", ".join(sorted(avail_suites))
                raise RuntimeError(f'unrecognized suite "{s}". Valid suites are: {valid_suites}')

    return list(suites)


@click.command()
@click.argument(
    "hictk-bin",
    type=click.Path(
        exists=True,
        executable=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=pathlib.Path,
    ),
    required=True,
)
@click.argument(
    "config-file",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        path_type=pathlib.Path,
    ),
    required=True,
)
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=pathlib.Path),
    help="Path to the folder with the test files.",
)
@click.option(
    "--suites",
    type=str,
    default="all",
    help="Comma-separated list of names of the tests to be executed.\n"
    "Should be one of:\n" + "\n - ".join(get_test_names(include_all=True)),
)
@click.option(
    "--threads",
    help="Specify the maximum number of CPU threads to be used.",
    type=click.IntRange(1, nproc()),
    default=1,
    show_default=True,
)
@click.option(
    "--verbosity",
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    default="info",
    help="Set verbosity level.",
)
@click.option(
    "--result-file",
    help="Path where to write the test results.",
    type=pathlib.Path,
)
@click.option(
    "--force",
    help="Force overwrite existing output file(s).",
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--no-cleanup",
    help="Don't clean temporary files upon exit.",
    default=False,
    is_flag=True,
    show_default=True,
)
def main(
    hictk_bin: pathlib.Path,
    data_dir: pathlib.Path,
    config_file: pathlib.Path,
    suites: str,
    threads: int,
    verbosity: str,
    result_file: pathlib.Path,
    force: bool,
    no_cleanup: bool,
):
    """
    Run hictk integration test suite.

    HICTK_BIN:   Path to hictk's binary.
    CONFIG_FILE: Path to the config.toml.
    """
    logging.basicConfig(level=parse_log_lvl(verbosity))

    if result_file and result_file.exists():
        if force:
            result_file.unlink()
        else:
            raise RuntimeError(f'refusing to ovrewrite file "{result_file}"')

    suites = parse_test_suites(suites)
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = init_results(hictk_bin)

    with WorkingDirectory(delete=not no_cleanup) as wd:
        hictk_bin = wd.stage_file(hictk_bin)
        for test in suites:
            mod = importlib.import_module(f"hictk_integration_suite.cli.{test}")

            t0 = time.time()
            config = import_config_and_stage_files(config_file, data_dir, wd, command=test)
            delta = (time.time() - t0) * 1000.0
            logging.info(f"staging test files for {test} tests took {delta:.2f}ms")

            t0 = time.time()
            plans = mod.plan_tests(hictk_bin, config, wd, threads)
            delta = (time.time() - t0) * 1000.0
            logging.info(f"planning for {test} tests took {delta:.2f}ms")

            t0 = time.time()
            res = mod.run_tests(plans, wd, no_cleanup)
            delta = time.time() - t0
            logging.info(f"running tests for {test} took {delta:.2f}s")

            num_pass += res[0]
            num_fail += res[1]
            num_skip += res[2]
            results["results"] |= res[3]

    results["results"]["pass"] = num_pass
    results["results"]["fail"] = num_fail
    results["results"]["fail"] = num_skip
    results["success"] = num_fail == 0
    if result_file is not None:
        with open(result_file, "w") as f:
            f.write(json.dumps(results, indent=2))

    print("", file=sys.stderr)
    print(f"# PASS: {num_pass}", file=sys.stderr)
    print(f"# SKIP: {num_skip}", file=sys.stderr)
    print(f"# FAIL: {num_fail}", file=sys.stderr)

    sys.exit(num_fail != 0)


if __name__ == "__main__":
    main()
