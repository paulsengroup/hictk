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
import tomllib
from typing import Any, Dict, List

import click

from hictk_integration_suite.runners.hictk.common import version


def get_test_names(include_all: bool = True) -> List[str]:
    if include_all:
        names = ["all"]
    else:
        names = []
    names += ["metadata"]
    return names


def update_uris(config: Dict, data_dir: pathlib.Path) -> Dict:
    def _update_uri(uri: pathlib.Path) -> str:
        path = os.path.join(data_dir, uri)
        if not os.path.exists(path):
            raise RuntimeError(f'file "{path}" does not exists')
        return path

    if "files" not in config:
        return config

    new_config = config.copy()
    for i, mappings in enumerate(config.get("files", [])):
        if "uri" in mappings:
            new_config["files"][i]["uri"] = _update_uri(mappings["uri"])

    return new_config


def import_config(path: pathlib.Path, data_dir: pathlib.Path | None, command: str | None = None) -> Dict[str, Any]:
    with open(path, "rb") as f:
        config = tomllib.load(f)

    if data_dir:
        for k, v in config.items():
            config[k] = update_uris(v, data_dir)

    if command is None:
        return config

    return config[command]


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
    raise NotImplemented()


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
    type=click.Choice(get_test_names(include_all=True)),
    default="all",
    help="Comma-separated list of names of the tests to be executed.\n"
    "Should be one of:\n" + "\n - ".join(get_test_names(include_all=True)),
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
)
def main(
    hictk_bin: pathlib.Path,
    data_dir: pathlib.Path,
    config_file: pathlib.Path,
    suites: str,
    verbosity: str,
    result_file: pathlib.Path,
):
    """
    Run hictk integration test suite.

    HICTK_BIN:   Path to hictk's binary.
    CONFIG_FILE: Path to the config.toml.
    """
    logging.basicConfig(level=parse_log_lvl(verbosity))

    if result_file and os.path.exists(result_file):
        raise RuntimeError(f'refusing to ovrewrite file "{result_file}"')

    suites = suites.split(",")
    if isinstance(suites, str):
        suites = [suites]
    else:
        suites = set(suites)

    if "all" in suites:
        suites = get_test_names(include_all=False)

    config = import_config(config_file, data_dir)

    num_pass = 0
    num_fail = 0
    results = init_results(hictk_bin)
    for test in suites:
        mod = importlib.import_module(f"hictk_integration_suite.cli.{test}")
        res = mod.run_tests(hictk_bin, config[test])
        num_pass += res[0]
        num_fail += res[1]
        results["results"] |= res[2]

    results["results"]["pass"] = num_pass
    results["results"]["fail"] = num_fail
    results["success"] = num_fail == 0
    if result_file is not None:
        with open(result_file, "w") as f:
            f.write(json.dumps(results, indent=2))

    print("")
    print(f"# PASS: {num_pass}")
    print(f"# FAIL: {num_fail}")
    sys.exit(num_fail != 0)


if __name__ == "__main__":
    main()
