# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
import shutil
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.common import parse_uri
from hictk_integration_suite.tests.rename_chromosomes import (
    HictkRenameChromosomes,
    HictkRenameChromosomesCli,
)

from .common import WorkingDirectory, _get_uri, _make_file_writeable, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    cool_uri: pathlib.Path,
    hic_uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-rename-chromosomes-cli",
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
        factory | {"args": tuple(("rename-chromosomes",))},
        factory | {"args": tuple(("rename-chromosomes", "--help")), "expect_failure": False},
        factory | {"args": tuple(("rename-chromosomes", "not-a-file"))},
        factory | {"args": tuple(("rename-chromosomes", "not-a-file", "--add-chr-prefix"))},
        factory | {"args": tuple(("rename-chromosomes", "--foobar"))},
        factory | {"args": tuple(("rename-chromosomes", str(cool_uri), "foobar"))},
        factory | {"args": tuple(("rename-chromosomes", str(cool_uri), "--foobar"))},
        factory | {"args": tuple(("rename-chromosomes", str(hic_uri), "--add-chr-prefix"))},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _stage_invalid_name_mappings(mappings: List[str], wd: WorkingDirectory) -> List[pathlib.Path]:
    paths = []
    parent_dir = wd.mkdtemp()
    for i, mapping in enumerate(mappings):
        path = wd.touch(parent_dir / f"invald_name_mappings_{i:03d}.txt")
        path.write_text(f"{mapping}\n")
        paths.append(path)

    return paths


def _stage_valid_name_mappings(
    mappings: List[Dict[str, str]], wd: WorkingDirectory
) -> Dict[pathlib.Path, pathlib.Path]:
    paths = {}
    parent_dir = wd.mkdtemp()
    for i, mapping in enumerate(mappings):
        payload = "\n".join(f"{chrom1}\t{chrom2}" for chrom1, chrom2 in mapping["mappings"].items())
        path = wd.touch(parent_dir / f"vald_name_mappings_{i:03d}.txt")
        path.write_text(f"{payload}\n")

        uri = wd[mapping["uri"]]
        paths[uri] = path

    return paths


def _stage_uri(uri: pathlib.Path | str, wd: WorkingDirectory) -> str:
    file_name, grp = parse_uri(uri)
    tmpdir = wd.mkdtemp()
    writeable_file = wd.touch(tmpdir / pathlib.Path(file_name).name)
    shutil.copy2(file_name, writeable_file)
    _make_file_writeable(writeable_file)

    if not grp:
        return str(writeable_file)

    return f"{writeable_file}::{grp}"


def _plan_tests_cmd(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-rename-chromosomes",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 1.0,
    }

    valid_name_mappings = _stage_valid_name_mappings(config["name-mappings"], wd)
    invalid_name_mappings = _stage_invalid_name_mappings(config["invalid-name-mappings"], wd)

    for c in config["test-cases"]:
        uri = wd[c["uri"]]
        expect_failure = c.get("expect-failure", False)
        for flag in ("--add-chr-prefix", "--remove-chr-prefix"):
            new_uri = _stage_uri(uri, wd)
            args = ["rename-chromosomes", str(new_uri), flag]
            plans.append(
                factory
                | {
                    "args": tuple(args),
                    "test_file": str(new_uri),
                    "expect_failure": expect_failure,
                }
            )

        new_uri = _stage_uri(uri, wd)
        plans.append(
            factory
            | {
                "args": (
                    "rename-chromosomes",
                    str(new_uri),
                    "--name-mappings",
                    str(valid_name_mappings[uri]),
                ),
                "test_file": str(new_uri),
                "name_mappings": str(valid_name_mappings[uri]),
                "expect_failure": expect_failure,
            }
        )

        for path in invalid_name_mappings:
            plans.append(
                factory
                | {
                    "args": (
                        "rename-chromosomes",
                        str(uri),
                        "--name-mappings",
                        str(path),
                    ),
                    "name_mappings": str(path),
                    "test_file": str(uri),
                    "expect_failure": True,
                }
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int = -1,
) -> List[ImmutableOrderedDict]:
    return _plan_tests_cli(hictk_bin, _get_uri(config, "cool"), _get_uri(config, "hic"), wd) + _plan_tests_cmd(
        hictk_bin, config, wd
    )


def run_tests(plans: List[ImmutableOrderedDict], wd: WorkingDirectory, no_cleanup: bool) -> Tuple[int, int, int, Dict]:
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
        assert title.startswith("hictk-rename-chromosomes")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkRenameChromosomesCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkRenameChromosomes(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    if not no_cleanup:
        wd.rmtree(cwd)
        wd.rmtree(tmpdir)

    return num_pass, num_fail, num_skip, results
