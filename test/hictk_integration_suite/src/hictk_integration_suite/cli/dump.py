# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from typing import Any, Dict, List, Set, Tuple

from hictk_integration_suite.tests.dump import HictkDump, HictkDumpCli


def _add_missing_reference_uris(config: Dict[str, Any]) -> Dict[str, Any]:
    for mapping in config["files"]:
        if "reference-uri" not in mapping:
            mapping["reference-uri"] = mapping["uri"]

    return config


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


def _merge_results(
    new_results: Tuple[int, int, Dict], old_results: Tuple[int, int, Dict] | None = None
) -> Tuple[int, int, Dict]:
    if old_results is None:
        return new_results

    old_results = list(old_results)
    old_results[0] += new_results[0]
    old_results[1] += new_results[1]
    old_results[2] |= new_results[2]

    return old_results[0], old_results[1], old_results[2]


def _test_help(hictk_bin: pathlib.Path, uri: pathlib.Path) -> Tuple[int, int, Dict]:
    test = HictkDumpCli(hictk_bin)

    num_pass = 0
    num_fail = 0
    results = []

    status = test.run(["dump"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["dump", "--help"], expect_failure=False)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["dump", "not-a-file"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["dump", "--foobar"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    status = test.run(["dump", str(uri), "--foobar"], expect_failure=True)
    results.append(status)
    num_pass += status["status"] == "PASS"
    num_fail += status["status"] == "FAIL"
    logging.info(status)

    return num_pass, num_fail, {"dump-cli": results}


def _test_invalid_args():
    pass


def _generate_test_args_for_hictk_dump_bins(config: Dict[str, Any]) -> Tuple:
    test_args = {}

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

            test_args |= {
                (tuple(query_gw), query["reference-uri"]): None,
                (tuple(query_subset1), query["reference-uri"]): None,
                (tuple(query_subset2), query["reference-uri"]): None,
            }

    logging.debug(f"hictk-dump-bins: generated {len(test_args)} test cases")

    return tuple(test_args.keys())


def _generate_test_args_for_hictk_dump_chroms(config: Dict[str, Any]) -> Tuple:
    test_args = {}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            args1 = _make_hictk_dump_args(
                query, drop_args={"range1", "range2", "normalization"}, add_args={"--table": "chroms"}
            )
            args2 = _make_hictk_dump_args(query, drop_args={"range2", "normalization"}, add_args={"--table": "chroms"})
            args3 = _make_hictk_dump_args(query, drop_args={"normalization"}, add_args={"--table": "chroms"})

            test_args |= {
                (tuple(args1), query["reference-uri"]): None,
                (tuple(args2), query["reference-uri"]): None,
                (tuple(args3), query["reference-uri"]): None,
            }

    logging.debug(f"hictk-dump-bins: generated {len(test_args)} test cases")

    return tuple(test_args.keys())


def _generate_test_args_for_hictk_dump_cis(config: Dict[str, Any]) -> Tuple:
    test_args = {}

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

            test_args |= {
                (tuple(query_raw_coo), query["reference-uri"]): None,
                (tuple(query_norm_coo), query["reference-uri"]): None,
                (tuple(query_cis_only_raw_coo), query["reference-uri"]): None,
                (tuple(query_cis_only_norm_coo), query["reference-uri"]): None,
                (tuple(query_raw_bg2), query["reference-uri"]): None,
                (tuple(query_norm_bg2), query["reference-uri"]): None,
                (tuple(query_cis_only_raw_bg2), query["reference-uri"]): None,
                (tuple(query_cis_only_norm_bg2), query["reference-uri"]): None,
            }

    logging.debug(f"hictk-dump-cis: generated {len(test_args)} test cases")

    return tuple(test_args.keys())


def _generate_test_args_for_hictk_dump_trans(config: Dict[str, Any]) -> Tuple:
    test_args = {}

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

            test_args |= {
                (tuple(query_raw_coo), query["reference-uri"]): None,
                (tuple(query_norm_coo), query["reference-uri"]): None,
                (tuple(query_cis_only_raw_coo), query["reference-uri"]): None,
                (tuple(query_cis_only_norm_coo), query["reference-uri"]): None,
                (tuple(query_raw_bg2), query["reference-uri"]): None,
                (tuple(query_norm_bg2), query["reference-uri"]): None,
                (tuple(query_trans_only_raw_bg2), query["reference-uri"]): None,
                (tuple(query_trans_only_norm_bg2), query["reference-uri"]): None,
            }

    logging.debug(f"hictk-dump-trans: generated {len(test_args)} test cases")
    return tuple(test_args.keys())


def _generate_test_args_for_hictk_dump_gw(config: Dict[str, Any]) -> Tuple:
    test_args = {}

    for c in config["files"]:
        for query in _extract_queries_for_uri(c["uri"], c["reference-uri"], c.get("resolution"), c.get("cell"), config):
            assert query.get("range1") is not None
            assert query.get("range2") is not None

            query_raw_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2", "normalization"})
            query_norm_coo = _make_hictk_dump_args(query, drop_args={"range1", "range2"})

            query_raw_bg2 = query_raw_coo + ["--join"]
            query_norm_bg2 = query_norm_coo + ["--join"]

            test_args |= {
                (tuple(query_raw_coo), query["reference-uri"]): None,
                (tuple(query_norm_coo), query["reference-uri"]): None,
                (tuple(query_raw_bg2), query["reference-uri"]): None,
                (tuple(query_norm_bg2), query["reference-uri"]): None,
            }

    logging.debug(f"hictk-dump-gw: generated {len(test_args)} test cases")
    return tuple(test_args.keys())


def _test_cmd(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> Tuple[int, int, Dict]:
    # We split test argument generation and test execution so that we can de-duplicate tests with identical
    # parameters, which can happen e.g. when the input file does not have balanced interactions
    test_cases = {
        "hictk-dump-bins": _generate_test_args_for_hictk_dump_bins(config),
        "hictk-dump-chroms": _generate_test_args_for_hictk_dump_chroms(config),
        "hictk-dump-pixels-cis": _generate_test_args_for_hictk_dump_cis(config),
        "hictk-dump-pixels-trans": _generate_test_args_for_hictk_dump_trans(config),
        "hictk-dump-pixels-gw": _generate_test_args_for_hictk_dump_gw(config),
    }

    timeouts = {
        "hictk-dump-bins": 5,
        "hictk-dump-chroms": 1,
        "hictk-dump-pixels-cis": 30,
        "hictk-dump-pixels-trans": 30,
        "hictk-dump-pixels-gw": 60,
    }

    test = HictkDump(hictk_bin)

    num_pass = 0
    num_fail = 0
    results = {}

    for title, test_args in test_cases.items():
        results[title] = []
        for args, reference in test_args:
            try:
                status = test.run(args, reference, timeout=timeouts[title], title=title)
                results[title].append(status)
                num_pass += status["status"] == "PASS"
                num_fail += status["status"] == "FAIL"
                logging.info(status)
            except:  # noqa
                logging.error(f"failed to executed test {[str(hictk_bin)] + list(args)}")
                raise

    return num_pass, num_fail, results


def run_tests(hictk_bin: pathlib.Path, config: Dict[str, Any]) -> Tuple[int, int, Dict]:
    config = _add_missing_reference_uris(config.copy())

    help_results = _test_help(hictk_bin, config["files"][0]["uri"])
    cmd_results = _test_cmd(hictk_bin, config)

    return _merge_results(help_results, cmd_results)
