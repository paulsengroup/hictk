# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, List, Set, Tuple

from immutabledict import ImmutableOrderedDict, immutabledict

from hictk_integration_suite.tests.dump import HictkDump, HictkDumpCli

from .common import _add_default_reference_uris, _get_uri, _preprocess_plan


def _extract_queries_for_uri(
    uri: str, reference_uri: str, resolution: int | None, cell: str | None, config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    queries = []
    for c in config["queries"]:
        if not uri.endswith(c["uri"]):
            continue

        if resolution:
            assert cell is None
            reference_uri = f"{reference_uri}::/resolutions/{resolution}"

        if cell:
            assert resolution is None
            uri = f"{uri}::/cells/{cell}"
            reference_uri = f"{reference_uri}::/cells/{cell}"

        queries.append(
            {
                "uri": uri,
                "reference-uri": reference_uri,
                "resolution": resolution,
                "range1": c.get("range1"),
                "range2": c.get("range2"),
                "normalization": c.get("normalization"),
            }
        )

    return queries


def _make_hictk_dump_args(
    config: Dict[str, Any], drop_args: Set[str] | None = None, add_args: Dict[str, Any] | None = None
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

    for k, v in add_args.items():
        if not v:
            args.append(str(k))
        else:
            args.extend((str(k), str(v)))

    return args


def _plan_tests_cli(
    hictk_bin: pathlib.Path, uri: pathlib.Path, title: str = "hictk-dump-cli"
) -> List[ImmutableOrderedDict]:
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0, "expect_failure": True}
    plans = (
        factory | {"args": tuple(("dump",))},
        factory | {"args": tuple(("dump", "--help")), "expect_failure": False},
        factory | {"args": tuple(("dump", "not-a-file"))},
        factory | {"args": tuple(("dump", str(uri), "foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--foobar"))},
        factory | {"args": tuple(("dump", str(uri), "--foobar"))},
    )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_bins(
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-dump-bins"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 5.0}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None

            query_gw = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--table": "bins"}
            )
            query_subset1 = _make_hictk_dump_args(
                query, drop_args={"range2", "normalization"}, add_args={"--table": "bins"}
            )
            query_subset2 = _make_hictk_dump_args(query, drop_args={"normalization"}, add_args={"--table": "bins"})

            plans.extend(
                (
                    factory | {"args": tuple(query_gw), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_subset1), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_subset2), "reference_uri": query["reference-uri"]},
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_chroms(
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-dump-chroms"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 1.0}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            args1 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--table": "chroms"}
            )
            args2 = _make_hictk_dump_args(query, drop_args={"range2", "normalization"}, add_args={"--table": "chroms"})
            args3 = _make_hictk_dump_args(query, drop_args={"normalization"}, add_args={"--table": "chroms"})

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
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-dump-cis"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 10.0}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None

            query_raw_coo = _make_hictk_dump_args(query, drop_args={"range2", "normalization"})
            query_norm_coo = _make_hictk_dump_args(query, drop_args={"range2"})
            query_cis_only_raw_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--cis-only": None}
            )
            query_cis_only_norm_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"--cis-only": None}
            )

            query_raw_bg2 = query_raw_coo + ["--join"]
            query_norm_bg2 = query_norm_coo + ["--join"]
            query_cis_only_raw_bg2 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--cis-only": None, "--join": None}
            )
            query_cis_only_norm_bg2 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"--cis-only": None, "--join": None}
            )

            plans.extend(
                (
                    factory | {"args": tuple(query_raw_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_raw_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_norm_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_raw_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_raw_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_norm_bg2), "reference_uri": query["reference-uri"]},
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_trans(
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-dump-trans"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 30.0}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            query_raw_coo = _make_hictk_dump_args(query, drop_args={"normalization"})
            query_norm_coo = _make_hictk_dump_args(query)
            query_cis_only_raw_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--trans-only": None}
            )
            query_cis_only_norm_coo = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"--trans-only": None}
            )

            query_raw_bg2 = query_raw_coo + ["--join"]
            query_norm_bg2 = query_norm_coo + ["--join"]
            query_trans_only_raw_bg2 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--trans-only": None}
            )
            query_trans_only_norm_bg2 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2"}, add_args={"--trans-only": None, "--join": None}
            )

            plans.extend(
                (
                    factory | {"args": tuple(query_raw_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_raw_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_cis_only_norm_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_raw_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_trans_only_raw_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_trans_only_norm_bg2), "reference_uri": query["reference-uri"]},
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests_hictk_dump_gw(
    hictk_bin: pathlib.Path, config: Dict[str, Any], title: str = "hictk-dump-gw"
) -> List[ImmutableOrderedDict]:
    plans = []
    factory = {"hictk_bin": str(hictk_bin), "title": title, "timeout": 45.0}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            query_raw_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2", "normalization"})
            query_norm_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2"})

            query_raw_bg2 = query_raw_coo + ["--join"]
            query_norm_bg2 = query_norm_coo + ["--join"]

            plans.extend(
                (
                    factory | {"args": tuple(query_raw_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_coo), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_raw_bg2), "reference_uri": query["reference-uri"]},
                    factory | {"args": tuple(query_norm_bg2), "reference_uri": query["reference-uri"]},
                )
            )

    plans = list(set(immutabledict(p) for p in plans))
    logging.debug(f"{title}: generated {len(plans)} test cases")
    return plans


def plan_tests(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> List[ImmutableOrderedDict]:
    config = _add_default_reference_uris(config.copy())
    return (
        _plan_tests_cli(hictk_bin, _get_uri(config))
        + _plan_tests_hictk_dump_bins(hictk_bin, config)
        + _plan_tests_hictk_dump_chroms(hictk_bin, config)
        + _plan_tests_hictk_dump_cis(hictk_bin, config)
        + _plan_tests_hictk_dump_trans(hictk_bin, config)
        + _plan_tests_hictk_dump_gw(hictk_bin, config)
    )


def run_tests(plans: List[ImmutableOrderedDict]) -> Tuple[int, int, int, Dict]:
    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    for p in plans:
        skip, p = _preprocess_plan(p)
        if skip:
            logging.info(f"SKIPPING {p}")
            num_skip += 1
            continue
        title = p["title"]
        assert title.startswith("hictk-dump")
        hictk = p.pop("hictk_bin")
        if title.endswith("-cli"):
            test = HictkDumpCli(hictk)
        else:
            test = HictkDump(hictk)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        logging.info(status)

    return num_pass, num_fail, num_skip, results
