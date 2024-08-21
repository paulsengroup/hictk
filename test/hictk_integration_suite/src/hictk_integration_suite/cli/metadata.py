# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.metadata import HictkMetadata, HictkMetadataCli

from .common import WorkingDirectory, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path, uri: pathlib.Path, wd: WorkingDirectory, title: str = "hictk-metadata-cli"
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0}
    plans = (
        factory | {"args": tuple(("metadata",)), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--help")), "expect_failure": False},
        factory | {"args": tuple(("metadata", "not-a-file")), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "--foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "--format", "foobar")), "expect_failure": True},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_cmd(
    hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory, title: str = "hictk-metadata"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0, "expect_failure": False}
    for c in config["files"]:
        factory["file_format"] = c["format"]
        factory["variable_bin_size"] = c.get("variable-bin-size", False)
        uri = wd[c["uri"]]
        for fmt in config["output-formats"]:
            for recursive in [True, False]:
                args = ["metadata", str(uri), "--output-format", fmt]
                if recursive:
                    args.append("--recursive")

                plans.append(factory | {"args": tuple(args)})

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config), wd) + _plan_tests_cmd(hictk_bin, config, wd)


def run_tests(plans: List[ImmutableOrderedDict], wd: WorkingDirectory) -> Tuple[int, int, int, Dict]:
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    cwd = wd.mkdtemp()
    tmpdir = wd.mkdtemp()

    for p in plans:
        skip, p = _preprocess_plan(p, wd)
        if skip:
            logging.info(f"SKIPPING {p}")
            num_skip += 1
            continue
        title = p["title"]
        assert title.startswith("hictk-metadata")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkMetadataCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkMetadata(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    return num_pass, num_fail, num_skip, results
