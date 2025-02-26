# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from typing import Any, Dict, List, Set, Tuple

from hictk_integration_suite.common import URI
from hictk_integration_suite.tests.dump import HictkDump, HictkDumpCli
from hictk_integration_suite.validators.file_formats import is_multires, is_scool
from immutabledict import ImmutableOrderedDict, immutabledict

from .common import WorkingDirectory, _argument_map_to_list, _get_uri, _preprocess_plan


def _extract_queries_for_uri(
    uri: pathlib.Path,
    reference_uri: pathlib.Path,
    resolution: int | None,
    cell: str | None,
    config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    file = str(URI(uri).path)
    queries = []
    for c in config["queries"]:
        if not file.endswith(os.path.basename(c["uri"])):
            continue

        queries.append(
            {
                "uri": str(uri),
                "reference-uri": str(reference_uri),
                "resolution": resolution,
                "range1": c.get("range1"),
                "range2": c.get("range2"),
                "normalization": c.get("normalization"),
            }
        )

    return queries


def _make_hictk_dump_args(
    config: Dict[str, Any],
    drop_args: Set[str] | None = None,
    add_args: Dict[str, Any] | None = None,
) -> List[str]:
    if drop_args is None:
        drop_args = set()
    if add_args is None:
        add_args = {}

    config = config.copy()
    for k in drop_args:
        config.pop(k)

    args = ["dump", config["uri"]]
    resolution = config.get("resolution")
    if resolution:
        args.extend(("--resolution", str(resolution)))
    range1 = config.get("range1")
    if range1:
        args.extend(("--range", range1))
    range2 = config.get("range2")
    if range2:
        args.extend(("--range2", range2))
    normalization = config.get("normalization")
    if normalization:
        args.extend(("--balance", normalization))
    join = config.get("join")
    if join:
        args.append("--join")

    args.extend(_argument_map_to_list(add_args))

    return args


def _plan_tests_cli(
    hictk_bin: pathlib.Path,
    uri: pathlib.Path,
    wd: WorkingDirectory,
    title: str = "hictk-dump-cli",
) -> List[ImmutableOrderedDict]:
    uri = wd[uri]
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
    }
    plans = (
        factory | {"args": tuple(("dump",))},
        factory | {"args": tuple(("dump", "--help")), "expect_failure": False},
        factory | {"args": tuple(("dump", "not-a-file"))},
        factory | {"args": tuple(("dump", str(uri), "foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--matrix-type", "foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--matrix-unit", "foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--table", "foobar"))},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_bins(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-dump-bins",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 5.0}

    for c in config["files"]:
        uri = wd[c["uri"]]
        reference_uri = wd[c.get("reference-uri", c["uri"])]
        resolution = c.get("resolution")
        cell = c.get("cell")
        factory["expect_failure"] = is_multires(uri) and resolution is None
        for query in _extract_queries_for_uri(uri, reference_uri, resolution, cell, config):
            assert query.get("range1") is not None

            # hictk dump ... -t bins
            query_gw = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"table": "bins"},
            )
            # hictk dump ... -t bins --range xxx
            query_subset1 = _make_hictk_dump_args(
                query, drop_args={"range2", "normalization"}, add_args={"table": "bins"}
            )
            # hictk dump ... -t bins --range xxx --range2 xxx
            query_subset2 = _make_hictk_dump_args(query, drop_args={"normalization"}, add_args={"table": "bins"})

            plans.extend(
                (
                    factory
                    | {
                        "args": tuple(query_gw),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_subset1),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_subset2),
                        "reference_uri": query["reference-uri"],
                    },
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_chroms(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-dump-chroms",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {
        "hictk_bin": str(hictk_bin),
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
    }

    for c in config["files"]:
        uri = wd[c["uri"]]
        reference_uri = wd[c.get("reference-uri", c["uri"])]
        for query in _extract_queries_for_uri(uri, reference_uri, c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            # hictk dump ... -t chroms
            args1 = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"table": "chroms"},
            )
            # hictk dump ... -t chroms --range xxx
            args2 = _make_hictk_dump_args(
                query,
                drop_args={"range2", "normalization"},
                add_args={"table": "chroms"},
            )
            # hictk dump ... -t chroms --range xxx --range2 xxx
            args3 = _make_hictk_dump_args(query, drop_args={"normalization"}, add_args={"table": "chroms"})

            plans.extend(
                (
                    factory | {"args": tuple(args1), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(args2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(args3), "reference_uri": query["reference-uri"]},
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_cis(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-dump-cis",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 10.0}

    for c in config["files"]:
        uri = wd[c["uri"]]
        reference_uri = wd[c.get("reference-uri", c["uri"])]
        factory["expect_failure"] = (is_multires(uri) and c.get("resolution") is None) or is_scool(uri)
        for query in _extract_queries_for_uri(uri, reference_uri, c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None

            # hictk dump ... --range xxx
            query_raw_coo = _make_hictk_dump_args(query, drop_args={"range2", "normalization"})
            # hictk dump ... --range xxx --normalization xxx
            query_norm_coo = _make_hictk_dump_args(query, drop_args={"range2"})

            # hictk dump ... --cis-only
            query_cis_only_raw_coo = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"cis-only": None},
            )

            # hictk dump ... --cis-only --normalization xxx
            query_cis_only_norm_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"cis-only": None}
            )

            # hictk dump ... --join
            query_raw_bg2 = query_raw_coo + ["--join"]
            # hictk dump ... --join --normalization xxx
            query_norm_bg2 = query_norm_coo + ["--join"]

            # hictk dump ... --cis-only --join
            query_cis_only_raw_bg2 = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"cis-only": None, "join": None},
            )
            # hictk dump ... --cis-only --join --normalization xxx
            query_cis_only_norm_bg2 = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2"},
                add_args={"cis-only": None, "join": None},
            )

            plans.extend(
                (
                    factory
                    | {
                        "args": tuple(query_raw_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_cis_only_raw_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_cis_only_norm_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_raw_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_cis_only_raw_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_cis_only_norm_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_trans(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-dump-trans",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 30.0}

    for c in config["files"]:
        uri = wd[c["uri"]]
        reference_uri = wd[c.get("reference-uri", c["uri"])]
        factory["expect_failure"] = (is_multires(uri) and c.get("resolution") is None) or is_scool(uri)
        for query in _extract_queries_for_uri(uri, reference_uri, c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            # hictk dump ... --range xxx --range2 xxx
            query_raw_coo = _make_hictk_dump_args(query, drop_args={"normalization"})
            # hictk dump ... --range xxx --range2 xxx --normalization xxx
            query_norm_coo = _make_hictk_dump_args(query)

            # hictk dump ... --trans-only
            query_trans_only_raw_coo = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"trans-only": None},
            )
            # hictk dump ... --trans-only --normalization xxx
            query_trans_only_norm_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"trans-only": None}
            )

            # hictk dump ... --range xxx --range2 xxx --join
            query_raw_bg2 = query_raw_coo + ["--join"]
            # hictk dump ... --range xxx --range2 xxx --join --normalization xxx
            query_norm_bg2 = query_norm_coo + ["--join"]

            # hictk dump ... --trans-only --join
            query_trans_only_raw_bg2 = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2", "normalization"},
                add_args={"trans-only": None},
            )
            # hictk dump ... --trans-only --join --normalization xxx
            query_trans_only_norm_bg2 = _make_hictk_dump_args(
                query,
                drop_args={"range1", "range2"},
                add_args={"trans-only": None, "join": None},
            )

            plans.extend(
                (
                    factory
                    | {
                        "args": tuple(query_raw_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_trans_only_raw_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_trans_only_norm_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_raw_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_trans_only_raw_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_trans_only_norm_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_gw(
    hictk_bin: pathlib.Path,
    config: Dict[str, Any],
    wd: WorkingDirectory,
    title: str = "hictk-dump-gw",
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 45.0}

    for c in config["files"]:
        uri = wd[c["uri"]]
        reference_uri = wd[c.get("reference-uri", c["uri"])]
        factory["expect_failure"] = (is_multires(uri) and c.get("resolution") is None) or is_scool(uri)
        for query in _extract_queries_for_uri(uri, reference_uri, c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            # hictk dump ...
            query_raw_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2", "normalization"})
            # hictk dump ... --normalization xxx
            query_norm_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2"})

            # hictk dump ... --join
            query_raw_bg2 = query_raw_coo + ["--join"]
            # hictk dump ... --join --normalization xxx
            query_norm_bg2 = query_norm_coo + ["--join"]

            plans.extend(
                (
                    factory
                    | {
                        "args": tuple(query_raw_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_coo),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_raw_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                    factory
                    | {
                        "args": tuple(query_norm_bg2),
                        "reference_uri": query["reference-uri"],
                    },
                )
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
    return (
        _plan_tests_cli(hictk_bin, _get_uri(config), wd)
        + _plan_tests_hictk_dump_bins(hictk_bin, config, wd)
        + _plan_tests_hictk_dump_chroms(hictk_bin, config, wd)
        + _plan_tests_hictk_dump_cis(hictk_bin, config, wd)
        + _plan_tests_hictk_dump_trans(hictk_bin, config, wd)
        + _plan_tests_hictk_dump_gw(hictk_bin, config, wd)
    )


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
        assert title.startswith("hictk-dump")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkDumpCli(hictk, cwd=cwd, tmpdir=tmpdir)
        else:
            test = HictkDump(hictk, cwd=cwd, tmpdir=tmpdir)

        status = test.run(**p, max_attempts=max_attempts)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    if not no_cleanup:
        wd.rmtree(cwd)
        wd.rmtree(tmpdir)

    return num_pass, num_fail, num_skip, results
