# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, Tuple

from hictk_integration_suite.tests.metadata import HictkMetadata, HictkMetadataCli


def _test_help(hictk_bin: pathlib.Path, uri: pathlib.Path) -> Tuple[int, int, Dict]:
    test = HictkMetadataCli(hictk_bin)

    num_pass = 0
    num_fail = 0
    results = []

    status = test.run(["metadata"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["metadata", "--help"], expect_failure=False)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["metadata", "not-a-file"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["metadata", "--foobar"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["metadata", str(uri), "--foobar"], expect_failure=True)
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

    for file in config["files"]:
        for fmt in config["output-formats"]:
            for recursive in [True, False]:
                args = ["metadata", file, "--output-format", fmt]
                if recursive:
                    args.append("--recursive")
                status = test.run(args, timeout=5)
                results.append(status)
                num_pass += status["status"] == "PASS"
                num_fail += status["status"] == "FAIL"
                logging.info(status)

    return num_pass, num_fail, {"metadata": results}


def run_tests(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> Tuple[int, int, Dict]:
    num_pass, num_fail, cli_results = _test_help(hictk_bin, config["files"][0])
    res = _test_cmd(hictk_bin, config)

    num_pass += res[0]
    num_fail += res[1]
    results = cli_results | res[2]

    return num_pass, num_fail, results
