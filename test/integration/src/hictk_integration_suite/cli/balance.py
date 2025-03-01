# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
import shutil
from typing import Any, Dict, List, Tuple

from hictk_integration_suite.common import URI
from hictk_integration_suite.tests.balance import (
    HictkBalanceICE,
    HictkBalanceICECli,
    HictkBalanceSCALE,
    HictkBalanceSCALECli,
    HictkBalanceVC,
    HictkBalanceVCCli,
)
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import (
    WorkingDirectory,
    _argument_map_to_list,
    _get_uri,
    _make_file_writeable,
    _preprocess_plan,
)


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    algorithm: str,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str | None = None,
) -> List[ImmutableOrderedDict]:
    if title is None:
        title = f"hictk-balance-{algorithm}-cli"

    uri = wd[uri]

    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("balance",))},
        factory | {"args": tuple(("balance", algorithm))},
        factory | {"args": tuple(("balance", "--help")), "expect_failure": False},
        factory | {"args": tuple(("balance", algorithm, "--help")), "expect_failure": False},
        factory | {"args": tuple(("balance", algorithm, "not-a-file"))},
        factory | {"args": tuple(("balance", algorithm, "--foobar"))},
        factory | {"args": tuple(("balance", algorithm, str(uri), "--foobar"))},
        factory | {"args": tuple(("balance", algorithm, str(uri), "--mode=invalid"))},
        factory | {"args": tuple(("balance", algorithm, str(uri), "--tmpdir", "not-a-folder"))},
        factory
        | {
            "args": tuple(("balance", str(uri))),
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
    algorithm: str,
    mode: str,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    threads: int,
    title: str | None = None,
) -> List[ImmutableOrderedDict]:
    assert threads > 0

    if title is None:
        title = f"hictk-balance-{algorithm}"

    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 120.0}
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = URI(c["uri"]).path
        reference = str(wd[c.get("reference-uri", wd[c["uri"]])])

        # Make a writeable copy of the input file
        dest = str(wd.touch(cwd / input_file.name))
        shutil.copy2(input_file, dest)
        input_file = dest
        _make_file_writeable(input_file)

        args = [
            "balance",
            algorithm,
            input_file,
            "--mode",
            mode,
            "--force",
        ]

        if algorithm != "vc":
            args.extend(
                (
                    "--tmpdir",
                    str(tmpdir),
                    "--compression-lvl",
                    "1",
                    "--chunk-size",
                    "100",
                    "--threads",
                    str(threads),
                )
            )

            args.extend(_argument_map_to_list(c.get("args", {})))

        plans.append(
            factory
            | {
                "args": tuple(args),
                "reference_file": reference,
                "test_file": input_file,
                "cwd": str(cwd),
                "expect_failure": c.get("expect-failure", False),
                "no_validate_weights": c.get("no-validate-weights", False),
            }
        )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(
    hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory, threads: int
) -> List[ImmutableOrderedDict]:
    plans = (
        _plan_tests_cli(hictk_bin, "ice", _get_uri(config, "cool"), wd)
        + _plan_tests_cli(hictk_bin, "scale", _get_uri(config, "cool"), wd)
        + _plan_tests_cli(hictk_bin, "vc", _get_uri(config, "cool"), wd)
    )

    for mode in ("gw", "cis", "trans"):
        plans.extend(_plan_tests_cmd(hictk_bin, "ice", mode, config, wd, threads))
        plans.extend(_plan_tests_cmd(hictk_bin, "scale", mode, config, wd, threads))
        plans.extend(_plan_tests_cmd(hictk_bin, "vc", mode, config, wd, threads))

    return plans


def _get_tester(
    hictk_bin: pathlib.Path, cwd: pathlib.Path, tmpdir: pathlib.Path, title: str
) -> (
    HictkBalanceICE | HictkBalanceSCALE | HictkBalanceVC | HictkBalanceICECli | HictkBalanceSCALECli | HictkBalanceVCCli
):
    if title.endswith("-ice-cli"):
        return HictkBalanceICECli(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-scale-cli"):
        return HictkBalanceSCALECli(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-vc-cli"):
        return HictkBalanceVCCli(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-ice"):
        return HictkBalanceICE(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-scale"):
        return HictkBalanceSCALE(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-vc"):
        return HictkBalanceVC(hictk_exec=hictk_bin, cwd=cwd, tmpdir=tmpdir)

    raise NotImplementedError


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
        assert title.startswith("hictk-balance")
        hictk = p.pop("hictk_bin")
        cwd = p.pop("cwd", main_cwd)
        test = _get_tester(hictk, title=title, cwd=cwd, tmpdir=tmpdir)
        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

        if not no_cleanup and cwd != main_cwd:
            wd.rmtree(cwd)

    if not no_cleanup:
        wd.rmtree(main_cwd)

    return num_pass, num_fail, num_skip, results
