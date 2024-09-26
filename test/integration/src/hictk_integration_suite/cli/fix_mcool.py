# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.fix_mcool import HictkFixMcool, HictkFixMcoolCli

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-fix-mcool-cli",
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0, "expect_failure": True}
    plans = (
        factory | {"args": tuple(("fix-mcool",))},
        factory | {"args": tuple(("fix-mcool", "--help")), "expect_failure": False},
        factory | {"args": tuple(("fix-mcool", "--foobar"))},
        factory | {"args": tuple(("fix-mcool", str(uri)))},
        factory | {"args": tuple(("fix-mcool", str(uri), "--foobar"))},
        factory | {"args": tuple(("fix-mcool", str(uri), "test.mcool", "--foobar"))},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_cmd(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int,
    title: str = "hictk-fix-mcool",
) -> List[ImmutableOrderedDict]:
    assert threads > 0
    plans = []
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "expect_failure": False,
    }
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = str(wd[c["input-uri"]])
        output_file = cwd / c["output"]
        expect_failure = c.get("expect-failure", False)
        timeout = c.get("timeout", 1.0)

        args = [
            "fix-mcool",
            input_file,
            str(output_file),
            "--tmpdir",
            str(tmpdir),
            "--threads",
            str(threads),
            "--compression-lvl",
            "1",
        ]
        args.extend(_argument_map_to_list(c.get("args", {})))
        plans.append(
            factory
            | {"args": tuple(args), "test_file": str(output_file), "expect_failure": expect_failure, "timeout": timeout}
        )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int,
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config, "mcool"), wd) + _plan_tests_cmd(hictk_bin, config, wd, threads)


def run_tests(
    plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool, max_attempts: int
) -> Tuple[int, int, int, Dict]:
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
        assert title.startswith("hictk-fix-mcool")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkFixMcoolCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkFixMcool(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    if not no_cleanup:
        wd.rmtree(cwd)
        wd.rmtree(tmpdir)

    return num_pass, num_fail, num_skip, results
