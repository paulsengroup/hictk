# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from typing import Any, Dict, List, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.balance import (
    HictkBalanceICE,
    HictkBalanceICECli,
    HictkBalanceSCALE,
    HictkBalanceSCALECli,
    HictkBalanceVC,
    HictkBalanceVCCli,
)

from .common import WorkingDirectory, _get_uri, _preprocess_plan


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    mode: str,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str | None = None,
) -> List[ImmutableOrderedDict]:
    if title is None:
        title = f"hictk-balance-{mode}-cli"

    uri = wd[uri]

    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 1.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("balance",))},
        factory | {"args": tuple(("balance", mode))},
        factory | {"args": tuple(("balance", "--help")), "expect_failure": False},
        factory | {"args": tuple(("balance", mode, "--help")), "expect_failure": False},
        factory | {"args": tuple(("balance", mode, "not-a-file"))},
        factory | {"args": tuple(("balance", mode, "--foobar"))},
        factory | {"args": tuple(("balance", mode, str(uri), "--foobar"))},
        factory | {"args": tuple(("balance", mode, str(uri), "--mode=invalid"))},
        factory | {"args": tuple(("balance", mode, str(uri), "--tmpdir", "not-a-folder"))},
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
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-balance",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 60.0}
    for c in config["test-cases"]:
        cwd = wd.mkdtemp(wd.name / title)
        tmpdir = wd.mkdir(cwd / "tmp")

        input_file = str(wd[c["input-uri"]])
        output_file = cwd / c["output"]
        reference = wd[c.get("reference-uri", input_file)]

        args = [
            "balance",
            str(input_file),
            str(output_file),
            "--tmpdir",
            str(tmpdir),
            "--compression-lvl",
            "1",
        ]
        for k, v in c.get("args", {}).items():
            k = "--" + (str(k).removeprefix("--"))
            if not v:
                args.append(k)
            elif isinstance(v, list):
                args.append(k)
                args.extend((str(x) for x in v))
            else:
                args.extend((k, str(v)))

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


def plan_tests(hictk_bin: pathlib.Path, config: Dict[str, Any], wd: WorkingDirectory) -> List[ImmutableOrderedDict]:
    return (
        _plan_tests_cli(hictk_bin, "ice", _get_uri(config, "cool"), wd)
        + _plan_tests_cli(hictk_bin, "scale", _get_uri(config, "cool"), wd)
        + _plan_tests_cli(hictk_bin, "vc", _get_uri(config, "cool"), wd)
    )


def _get_tester(
    hictk_bin: pathlib.Path, cwd: pathlib.Path, tmpdir: pathlib.Path, title: str
) -> (
    HictkBalanceICE | HictkBalanceSCALE | HictkBalanceVC | HictkBalanceICECli | HictkBalanceSCALECli | HictkBalanceVCCli
):
    if title.endswith("-ice"):
        return HictkBalanceICE(hictk_exec=hictk_bin, algorithm="ice", cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-scale"):
        return HictkBalanceSCALE(hictk_exec=hictk_bin, algorithm="scale", cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-vc"):
        return HictkBalanceVC(hictk_exec=hictk_bin, algorithm="vc", cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-ice-cli"):
        return HictkBalanceICECli(hictk_exec=hictk_bin, algorithm="ice", cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-scale-cli"):
        return HictkBalanceSCALECli(hictk_exec=hictk_bin, algorithm="scale", cwd=cwd, tmpdir=tmpdir)
    if title.endswith("-vc-cli"):
        return HictkBalanceVCCli(hictk_exec=hictk_bin, algorithm="vc", cwd=cwd, tmpdir=tmpdir)

    raise NotImplementedError


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
        assert title.startswith("hictk-balance")
        hictk = p.pop("hictk_bin")
        cwd = p.pop("cwd", main_cwd)
        test = _get_tester(hictk, title=title, cwd=cwd, tmpdir=tmpdir)
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
