# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import pathlib
import platform
from typing import Any, Dict, List, Mapping, Tuple


def _add_default_reference_uris(config: Mapping[str, Any]) -> Dict[str, Any]:
    config = dict(config)
    for mapping in config["files"]:
        if "reference-uri" not in mapping:
            mapping["reference-uri"] = mapping["uri"]

    return config


def _check_if_test_should_run(config: Mapping[str, Any]) -> bool:
    if config.get("skip_windows") and platform.system() == "Windows":
        return False
    if config.get("skip_linux") and platform.system() == "Linux":
        return False
    if config.get("skip_macos") and platform.system() == "Darwin":
        return False
    if config.get("skip_unix") and platform.system() in {"Linux", "Darwin"}:
        return False

    return True


def _strip_fields_from_config(config: Mapping[str, Any], fields: List[str] | None = None) -> Dict[str, Any]:
    if fields is None:
        fields = [
            "skip_windows",
            "skip_linux",
            "skip_macos",
            "skip_unix",
        ]

    config = dict(config)
    for f in fields:
        config.pop(f, None)

    return config


def _preprocess_plan(plan: Mapping[str, Any], fields_to_strip: List[str] | None = None) -> Tuple[bool, Dict[str, Any]]:
    return not _check_if_test_should_run(plan), _strip_fields_from_config(plan, fields_to_strip)


def _get_uri(config: Dict[str, Any], fmt: str | None = None) -> pathlib.Path:
    if "files" not in config or len(config["files"]) == 0:
        raise ValueError("unable to fetch uri from config")

    if fmt is None:
        return pathlib.Path(config["files"][0]["uri"])

    for c in config["files"]:
        if c["format"] == fmt:
            return pathlib.Path(c["uri"])

    raise ValueError(f'unable to fetch uri with format "{fmt}" from config')
