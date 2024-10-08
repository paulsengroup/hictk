# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from typing import Any, Dict, List, Tuple

from hictk_integration_suite.tests.load import HictkLoad, HictkLoadCli
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    pairs: pathlib.Path,
    chrom_sizes: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-load-cli",
) -> List[ImmutableOrderedDict]:
    pairs = str(wd[pairs])
    chrom_sizes = str(wd[chrom_sizes])
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("load",))},
        factory | {"args": tuple(("load", "--help")), "expect_failure": False},
        factory | {"args": tuple(("load", pairs))},
        factory | {"args": tuple(("load", "not-a-file1", "--chrom-sizes", "not-a-file2", "--bin-size", "10000"))},
        factory
        | {"args": tuple(("load", pairs, "not-a-file1", "--chrom-sizes", "not-a-file2", "--bin-size", "10000"))},
        factory | {"args": tuple(("load", pairs, "--foobar"))},
        factory
        | {
            "args": tuple(
                (
                    "load",
                    pairs,
                    "test.cool",
                    "--chrom-sizes",
                    chrom_sizes,
                    "--bin-size",
                    "10000",
                    "--tmpdir",
                    "not-a-folder",
                )
            )
        },
        factory
        | {
            "args": tuple(("load", pairs, "test.cool", "--chrom-sizes", chrom_sizes, "--bin-size", "10000")),
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
    title: str = "hictk-load",
) -> List[ImmutableOrderedDict]:
    assert threads > 0
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 60.0}
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = str(wd[c["input-path"]])

        chrom_sizes = wd.get(c.get("chrom-sizes-path"))
        bins = wd.get(c.get("bin-table-path"))

        for output_name in config["output-names"]:
            expect_failure = c.get("expect-failure", False)
            output_file = cwd / output_name
            args = [
                "load",
                input_file,
                str(output_file),
                "--tmpdir",
                str(tmpdir),
                "--compression-lvl",
                "1",
                "--threads",
                str(threads),
                "--chunk-size",
                "250000",
            ]

            if chrom_sizes is not None:
                assert bins is None
                args.extend(("--chrom-sizes", str(chrom_sizes)))

            if bins is not None:
                assert chrom_sizes is None
                args.extend(("--bin-table", str(bins)))
                if output_file.suffix == ".hic":
                    expect_failure = True

            args.extend(_argument_map_to_list(c.get("args", {})))

            reference_uri = wd.get(c.get("reference-uri"))
            if reference_uri is not None:
                reference_uri = str(reference_uri)

            try:
                i = args.index("--bin-size")
                resolution = int(args[i + 1])
            except ValueError:
                resolution = None

            plans.append(
                factory
                | {
                    "args": tuple(args),
                    "test_file": str(output_file),
                    "reference_uri": reference_uri,
                    "resolution": resolution,
                    "no_validate": c.get("no-validate", False),
                    "cwd": str(cwd),
                    "expect_failure": expect_failure,
                    "timeout": c.get("timeout", 90.0),
                }
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def get_chrom_sizes(config: Dict[str, Any]) -> pathlib.Path:
    if "files" not in config or len(config["files"]) == 0:
        raise ValueError("unable to fetch uri from config")

    for c in config["files"]:
        if c["format"] == "chrom.sizes":
            return pathlib.Path(c["input-path"])

    raise ValueError(f'unable to fetch uri with format "chrom.sizes" from config')


def get_pairs(config: Dict[str, Any]) -> pathlib.Path:
    if "files" not in config or len(config["files"]) == 0:
        raise ValueError("unable to fetch uri from config")

    for c in config["files"]:
        if c["format"] == "pairs":
            return pathlib.Path(c["input-path"])

    raise ValueError(f'unable to fetch uri with format "pairs" from config')


def plan_tests(
    hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory, threads: int
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, get_pairs(config), get_chrom_sizes(config), wd) + _plan_tests_cmd(
        hictk_bin, config, wd, threads
    )


def run_tests(
    plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool, max_attempts: int
) -> Tuple[int, int, int, Dict]:
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
        assert title.startswith("hictk-load")
        hictk = p.pop("hictk_bin")
        cwd = p.pop("cwd", main_cwd)
        if title.endswith("-cli"):
            test = HictkLoadCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkLoad(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    if not no_cleanup:
        wd.rmtree(main_cwd)

    return num_pass, num_fail, num_skip, results
