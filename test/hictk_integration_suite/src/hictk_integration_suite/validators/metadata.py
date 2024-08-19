# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import json as _json
import pathlib
import tomllib as _toml
from typing import Any, Dict, List, Set

import yaml as _yaml


def _detect_missing_fields(data: Dict[str, Any], required_fields: Set[str]) -> List[str]:
    missing_fields = []
    for field in required_fields:
        if field not in data:
            missing_fields.append(field)

    return missing_fields


def _parse_and_validate_json(handle_or_str, required_fields: Set[str]) -> bool:
    if isinstance(handle_or_str, str):
        data = _json.loads(handle_or_str)
    else:
        data = _json.load(handle_or_str)

    missing_fields = _detect_missing_fields(data, required_fields)

    if len(missing_fields) != 0:
        raise RuntimeError(f"JSON string is missing the following fields: {missing_fields}")

    return True


def _parse_and_validate_toml(handle_or_str, required_fields: Set[str]) -> bool:
    if isinstance(handle_or_str, str):
        data = _toml.loads(handle_or_str)
    else:
        data = _toml.load(handle_or_str)

    missing_fields = _detect_missing_fields(data, required_fields)

    if len(missing_fields) != 0:
        raise RuntimeError(f"TOML string is missing the following fields: {missing_fields}")

    return True


def _parse_and_validate_yaml(handle_or_str, required_fields: Set[str]) -> bool:
    data = _yaml.load(handle_or_str, _yaml.Loader)

    missing_fields = _detect_missing_fields(data, required_fields)

    if len(missing_fields) != 0:
        raise RuntimeError(f"YAML string is missing the following fields: {missing_fields}")

    return True


def json(
    handle_or_str,
    required_fields: Set[str] | None = None,
    raise_on_failure: bool = True,
) -> bool:
    if required_fields is None:
        required_fields = {}

    try:
        if isinstance(handle_or_str, pathlib.Path):
            with open(handle_or_str) as f:
                return _parse_and_validate_json(f, required_fields)
        else:
            return _parse_and_validate_json(handle_or_str, required_fields)
    except ValueError:
        if raise_on_failure:
            raise
        return False
    except RuntimeError:
        if raise_on_failure:
            raise
        return False


def toml(
    handle_or_str,
    required_fields: Set[str] | None = None,
    raise_on_failure: bool = True,
) -> bool:
    if required_fields is None:
        required_fields = {}

    try:
        if isinstance(handle_or_str, pathlib.Path):
            with open(handle_or_str) as f:
                return _parse_and_validate_toml(f, required_fields)
        else:
            return _parse_and_validate_toml(handle_or_str, required_fields)
    except ValueError:
        if raise_on_failure:
            raise
        return False
    except RuntimeError:
        if raise_on_failure:
            raise
        return False


def yaml(
    handle_or_str,
    required_fields: Set[str] | None = None,
    raise_on_failure: bool = True,
) -> bool:
    if required_fields is None:
        required_fields = {}

    try:
        if isinstance(handle_or_str, pathlib.Path):
            with open(handle_or_str) as f:
                return _parse_and_validate_yaml(f, required_fields)
        else:
            return _parse_and_validate_yaml(handle_or_str, required_fields)
    except ValueError:
        if raise_on_failure:
            raise
        return False
    except RuntimeError:
        if raise_on_failure:
            raise
        return False
