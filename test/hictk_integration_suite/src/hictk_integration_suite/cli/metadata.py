# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.metadata import HictkMetadata, HictkMetadataCli


def _plan_test_help(
    hictk_bin: pathlib.Path, uri: pathlib.Path, title: str = "hictk-metadata-cli"
) -> List[ImmutableOrderedDict]:
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0}
    plans = (
        factory | {"args": tuple(("metadata",)), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--help")), "expect_failure": False},
        factory | {"args": tuple(("metadata", "not-a-file")), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "--foobar")), "expect_failure": True},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_test_cmd(
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-metadata"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0, "expect_failure": False}
    for c in config["files"]:
        uri = c["uri"]
        for fmt in config["output-formats"]:
            for recursive in [True, False]:
                args = ["metadata", uri, "--output-format", fmt]
                if recursive:
                    args.append("--recursive")

                plans.append(factory | {"args": tuple(args)})

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _test_help(test_plans: List[Dict]) -> Tuple[int, int, Dict]:
    num_pass = 0
    num_fail = 0
    results = []

    for plan in test_plans:
        hictk_bin = plan.pop("hictk_bin")
        status = HictkMetadataCli(hictk_bin).run(**plan)
        results.append(status)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        logging.info(status)

    return num_pass, num_fail, {"metadata-cli": results}


def _test_cmd(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> Tuple[int, int, Dict]:
    test = HictkMetadata(hictk_bin)

    num_pass = 0
    num_fail = 0
    results = []

    for c in config["files"]:
        uri = c["uri"]
        for fmt in config["output-formats"]:
            for recursive in [True, False]:
                args = ["metadata", uri, "--output-format", fmt]
                if recursive:
                    args.append("--recursive")
                status = test.run(args, timeout=5)
                results.append(status)
                num_pass += status["status"] == "PASS"
                num_fail += status["status"] == "FAIL"
                logging.info(status)

    return num_pass, num_fail, {"metadata": results}


def plan_tests(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> List[ImmutableOrderedDict]:
    return _plan_test_help(hictk_bin, config["files"][0]["uri"]) + _plan_test_cmd(hictk_bin, config)


def run_tests(plans: List[ImmutableOrderedDict]) -> Tuple[int, int, Dict]:
    num_pass = 0
    num_fail = 0
    results = {}

    for p in plans:
        p = dict(p)
        title = p["title"]
        assert title.startswith("hictk-metadata")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkMetadataCli(hictk)
        else:
            test = HictkMetadata(hictk)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    return num_pass, num_fail, results
