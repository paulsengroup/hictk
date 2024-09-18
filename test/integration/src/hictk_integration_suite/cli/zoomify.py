# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.zoomify import HictkZoomify, HictkZoomifyCli

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-zoomify-cli",
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 1.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("zoomify",))},
        factory | {"args": tuple(("zoomify", "--help")), "expect_failure": False},
        factory | {"args": tuple(("zoomify", "not-a-file", "test.mcool"))},
        factory | {"args": tuple(("zoomify", "--foobar"))},
        factory | {"args": tuple(("zoomify", uri.as_posix(), "test.mcool", "foobar"))},
        factory | {"args": tuple(("zoomify", uri.as_posix(), "test.mcool", "--foobar"))},
        factory | {"args": tuple(("zoomify", uri.as_posix(), "test.mcool", "--tmpdir", "not-a-folder"))},
        factory
        | {
            "args": tuple(("zoomify", uri.as_posix(), "test.mcool")),
            "env_variables": immutabledict(
                os.environ | {var: "not-a-folder" for var in ("TMPDIR", "TMP", "TEMP", "TEMPDIR")}
            ),
            "skip_windows": True,
        },
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_cmd(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int,
    title: str = "hictk-zoomify",
) -> List[ImmutableOrderedDict]:
    assert threads > 0
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 60.0}
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = wd[c["input-uri"]].as_posix()
        output_file = cwd / c["output"]
        reference = wd[c.get("reference-uri", input_file)].as_posix()
        resolutions = c["resolutions"]

        args = [
            "zoomify",
            input_file,
            str(output_file),
            "--tmpdir",
            str(tmpdir),
            "--compression-lvl",
            "1",
            "--threads",
            str(threads),
            "--resolutions",
        ]

        args.extend(str(res) for res in resolutions)
        args.extend(_argument_map_to_list(c.get("args", {})))

        plans.append(
            factory
            | {
                "args": tuple(args),
                "reference_file": reference,
                "test_file": str(output_file),
                "resolutions": tuple(resolutions),
                "cwd": str(cwd),
                "expect_single_resolution": c.get("single-resolution", False),
                "expect_failure": c.get("expect-failure", False),
            }
        )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory, threads: int
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config), wd) + _plan_tests_cmd(hictk_bin, config, wd, threads)


def run_tests(plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool) -> Tuple[int, int, int, Dict]:
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    main_cwd = wd.mkdtemp()
    tmpdir = wd.mkdir(main_cwd / "tmp")

    for p in plans:
        skip, p = _preprocess_plan(p, wd)
        if skip:
            logging.info(f"SKIPPING {p}")
            num_skip += 1
            continue
        title = p["title"]
        assert title.startswith("hictk-zoomify")
        hictk = p.pop("hictk_bin")
        cwd = p.pop("cwd", main_cwd)
        if title.endswith("-cli"):
            test = HictkZoomifyCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkZoomify(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    if not no_cleanup:
        wd.rmtree(main_cwd)

    return num_pass, num_fail, num_skip, results
