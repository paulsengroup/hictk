# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Any, Dict, List, Tuple

import structlog
from hictk_integration_suite.tests.validate import HictkValidate, HictkValidateCli
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-validate-cli",
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 5.0, "expect_failure": True}
    plans = (
        factory | {"args": tuple(("validate",))},
        factory | {"args": tuple(("validate", "--help")), "expect_failure": False},
        factory | {"args": tuple(("validate", "--foobar"))},
        factory | {"args": tuple(("validate", str(uri), "foobar"))},
        factory | {"args": tuple(("validate", str(uri), "--foobar"))},
        factory | {"args": tuple(("validate", str(uri), "--format", "foobar"))},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_cmd(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-validate",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "expect_failure": False,
    }
    for c in config["test-cases"]:
        uri = str(wd[c["uri"]])
        expect_failure = c.get("expect-failure", False)
        timeout = c.get("timeout", 5.0)

        for fmt in config["output-formats"]:
            args = ["validate", uri, "--output-format", fmt]
            args.extend(_argument_map_to_list(c.get("args", {})))
            plans.append(factory | {"args": tuple(args), "expect_failure": expect_failure, "timeout": timeout})

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int = -1,
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config), wd) + _plan_tests_cmd(hictk_bin, config, wd)


def run_tests(
    plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool, max_attempts: int
) -> Tuple[int, int, int, Dict]:
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    cwd = wd.mkdtemp()
    tmpdir = wd.mkdtemp()

    logger = structlog.get_logger().bind()

    for p in plans:
        skip, p = _preprocess_plan(p, wd)
        if skip:
            logger.bind(status="SKIP").info(str(p))
            num_skip += 1
            continue
        title = p["title"]
        assert title.startswith("hictk-validate")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkValidateCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkValidate(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        if status["status"] == "PASS":
            logger.bind(**status).info("")
        else:
            logger.bind(**status).warning("")

    if not no_cleanup:
        wd.rmtree(cwd)
        wd.rmtree(tmpdir)

    return num_pass, num_fail, num_skip, results
