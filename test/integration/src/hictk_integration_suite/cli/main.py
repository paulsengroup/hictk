# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
from typing import Any, Dict, List, Tuple

import structlog
from hictk_integration_suite.tests.cli import HictkCli
from hictk_integration_suite.tests.verbosity import HictkVerbosity
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _get_uri, _preprocess_plan


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


def _plan_tests_verbosity(
    hictk_bin: pathlib.Path,
    uri: str,
    title: str = "hictk-verbosity",
) -> List[ImmutableOrderedDict]:
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
        "args": tuple(("dump", uri, "-t", "chroms")),
    }

    plans = (
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_QUIET": "1"})},
        factory | {"env_variables": immutabledict(os.environ | {"VERBOSE": "1"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "0"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "1"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "2"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "3"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "4"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "5"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "CRITICAL"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "ERROR"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "ERR"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "WARNING"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "WARN"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "INFO"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "DEBUG"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "critical"})},
        factory | {"env_variables": immutabledict(os.environ | {"HICTK_VERBOSITY": "foobar"})},
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
    return _plan_tests_cli(hictk_bin) + _plan_tests_verbosity(hictk_bin, _get_uri(config))


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
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkCli(hictk, cwd=main_cwd)
        else:
            test = HictkVerbosity(hictk, cwd=main_cwd)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        if status["status"] == "PASS":
            logger.bind(**status).info("")
        else:
            logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
