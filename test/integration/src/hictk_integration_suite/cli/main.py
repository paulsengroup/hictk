# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Any, Dict, List, Tuple

import structlog
from hictk_integration_suite.tests.cli import HictkCli
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    title: str = "hictk-cli",
) -> List[ImmutableOrderedDict]:
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
    }
    plans = (
        factory | {"args": tuple(("--help",))},
        factory | {"args": tuple(("--help-cite",))},
        factory | {"args": tuple(("--help-docs",))},
        factory | {"args": tuple(("--help-license",))},
        factory | {"args": tuple(("--help-telemetry",))},
        factory | {"args": tuple(("--foobar",)), "expect_failure": True},
        factory | {"args": tuple(("--help-cite", "--help-telemetry")), "expect_failure": True},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int = -1,
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin)


def run_tests(
    plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool, max_attempts: int
) -> Tuple[int, int, int, Dict]:
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    main_cwd = wd.mkdtemp()

    logger = structlog.get_logger().bind()

    for p in plans:
        skip, p = _preprocess_plan(p, wd)
        if skip:
            logger.bind(status="SKIP").info(str(p))
            num_skip += 1
            continue
        title = p["title"]
        assert title.startswith("hictk-cli")
        hictk = p.pop("hictk_bin")
        test = HictkCli(hictk, cwd=main_cwd)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        if status["status"] == "PASS":
            logger.bind(**status).info("")
        else:
            logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
