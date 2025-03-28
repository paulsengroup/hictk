# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Any, Dict, List, Tuple

import structlog
from hictk_integration_suite.tests.metadata import HictkMetadata, HictkMetadataCli
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-metadata-cli",
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 5.0}
    plans = (
        factory | {"args": tuple(("metadata",)), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--help")), "expect_failure": False},
        factory | {"args": tuple(("metadata", "not-a-file")), "expect_failure": True},
        factory | {"args": tuple(("metadata", "--foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "foobar")), "expect_failure": True},
        factory | {"args": tuple(("metadata", str(uri), "--foobar")), "expect_failure": True},
        factory
        | {
            "args": tuple(("metadata", str(uri), "--format", "foobar")),
            "expect_failure": True,
        },
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_cmd(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-metadata",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
    }
    for c in config["files"]:
        factory["file_format"] = c["format"]
        factory["variable_bin_size"] = c.get("variable-bin-size", False)
        uri = str(wd[c["uri"]])
        for fmt in config["output-formats"]:
            for recursive in [True, False]:
                args = ["metadata", uri, "--output-format", fmt]
                if recursive:
                    args.append("--recursive")

                plans.append(factory | {"args": tuple(args)})

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
        assert title.startswith("hictk-metadata")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkMetadataCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkMetadata(hictk, cwd=cwd, tmpdir=tmpdir)

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
