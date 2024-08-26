# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.convert import HictkConvert, HictkConvertCli

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    cool_uri: pathlib.Path,
    hic_uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-convert-cli",
) -> List[ImmutableOrderedDict]:
    cool_uri = wd[cool_uri]
    hic_uri = wd[hic_uri]

    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 1.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("convert",))},
        factory | {"args": tuple(("convert", "--help")), "expect_failure": False},
        factory | {"args": tuple(("convert", "not-a-file"))},
        factory | {"args": tuple(("convert", "not-a-file1", "not-a-file2"))},
        factory | {"args": tuple(("convert", "--foobar"))},
        factory | {"args": tuple(("convert", str(cool_uri), "foobar"))},
        factory | {"args": tuple(("convert", str(cool_uri), "--foobar"))},
        factory | {"args": tuple(("convert", str(cool_uri), "test.scool"))},
        factory | {"args": tuple(("convert", str(hic_uri), "test.scool"))},
        factory | {"args": tuple(("convert", str(cool_uri), "test.hic", "--tmpdir", "not-a-folder"))},
        factory | {"args": tuple(("convert", str(hic_uri), "test.mcool", "--tmpdir", "not-a-folder"))},
        factory
        | {
            "args": tuple(("convert", str(cool_uri), "test.hic")),
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
    title: str = "hictk-convert",
) -> List[ImmutableOrderedDict]:
    assert threads > 0
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 60.0}
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = str(wd[c["input-uri"]])
        output_file = cwd / c["output"]
        reference = wd[c.get("reference-uri", input_file)]

        args = [
            "convert",
            str(input_file),
            str(output_file),
            "--tmpdir",
            str(tmpdir),
            "--compression-lvl",
            "1",
            "--threads",
            str(threads),
        ]

        args.extend(_argument_map_to_list(c.get("args", {})))

        resolutions = c.get("args", {}).get("resolutions", [])
        plans.append(
            factory
            | {
                "args": tuple(args),
                "reference_file": str(reference),
                "test_file": str(output_file),
                "resolutions": tuple(resolutions),
                "cwd": str(cwd),
                "expect_failure": c.get("expect-failure", False),
            }
        )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory, threads: int
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config, "cool"), _get_uri(config, "hic"), wd) + _plan_tests_cmd(
        hictk_bin, config, wd, threads
    )


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
        assert title.startswith("hictk-convert")
        hictk = p.pop("hictk_bin")
        cwd = p.pop("cwd", main_cwd)
        if title.endswith("-cli"):
            test = HictkConvertCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkConvert(hictk, cwd=cwd, tmpdir=tmpdir)
        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

        if not no_cleanup and cwd != main_cwd:
            wd.rmtree(cwd)

    if not no_cleanup:
        wd.rmtree(main_cwd)

    return num_pass, num_fail, num_skip, results
