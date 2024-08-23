# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT
import hashlib
import json
import logging
import os.path
import pathlib
import platform
import shutil
import stat
import sys
import tempfile
from typing import Any, Dict, List, Mapping, Tuple

from immutabledict import immutabledict


class WorkingDirectory:
    def __init__(self, path: pathlib.Path | str | None = None, delete: bool = True):
        if path is None:
            path = tempfile.TemporaryDirectory(prefix="hictk-integration-test-", delete=False).name
        else:
            if pathlib.Path(path).exists():
                raise RuntimeError(f'"{path}" already exists')
            os.mkdir(path)

        self._delete = delete
        self._path = pathlib.Path(path)
        self._mappings = {}

    def __str__(self) -> str:
        return str(self._path)

    def __repr__(self) -> str:
        return f'WorkingDirectory("{self}")'

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    @staticmethod
    def _make_read_only(path: pathlib.Path | str):
        mode = stat.S_IRUSR
        if os.access(path, os.X_OK):
            mode |= stat.S_IXUSR
        pathlib.Path(path).chmod(mode)

    @staticmethod
    def _make_writeable(path: pathlib.Path | str):
        mode = stat.S_IRUSR | stat.S_IWUSR
        if os.access(path, os.X_OK):
            mode |= stat.S_IXUSR
        pathlib.Path(path).chmod(mode)

    @staticmethod
    def _parse_uri(s: pathlib.Path | str) -> Tuple[pathlib.Path, pathlib.Path | None]:
        path, _, grp = str(s).partition("::")
        if grp:
            grp = pathlib.Path(grp)
        else:
            grp = None

        return pathlib.Path(path), grp

    def stage_file(
        self,
        src: pathlib.Path | str,
        make_read_only: bool = True,
        exists_ok: bool = False,
    ) -> pathlib.Path:
        src = pathlib.Path(src).resolve()
        if src in self._mappings:
            return self._mappings[src]

        path, grp = self._parse_uri(src)
        if not path.exists():
            raise RuntimeError(f'source file "{path}" does not exist')

        dest_dir = self._path / "staged_files"
        dest_dir.mkdir(exist_ok=True)

        dest = dest_dir / path.name
        if dest.exists():
            if not exists_ok:
                raise RuntimeError(f'refusing to overwrite file "{dest}"')

            logging.debug(f'file "{path}" was already staged')

            if src != path:
                if grp:
                    uri = pathlib.Path(f"{self._mappings[path]}::{grp}")
                    self._mappings[src] = uri
                    self._mappings[uri] = uri
                else:
                    self._mappings[src] = self._mappings[path]
                    self._mappings[path] = self._mappings[path]

            return self._mappings[src]

        logging.debug(f'staging file "{path}"...')
        shutil.copy2(path, dest)
        if make_read_only:
            self._make_read_only(dest)

        dest = dest.resolve()
        self._mappings[path] = dest
        self._mappings[dest] = dest
        if src != path:
            if not grp:
                self._mappings[src] = dest
                self._mappings[dest] = dest
            else:
                uri = pathlib.Path(f"{dest}::{grp}")
                self._mappings[src] = uri
                self._mappings[uri] = uri

        return self._mappings[src]

    def mkdtemp(self, prefix: pathlib.Path | None = None) -> pathlib.Path:
        if prefix is None:
            prefix = self._path
        elif not self._path_belongs_to_wd(prefix, check_if_exists=False):
            raise RuntimeError(f'prefix "{prefix}" does not live under {self._path}')

        if not prefix.exists():
            prefix.mkdir()

        return pathlib.Path(tempfile.mkdtemp(dir=prefix))

    def tmpdir(self) -> pathlib.Path:
        path = self._path / "tmp"
        if not path.exists():
            path.mkdir()
        return path

    def mkdir(self, path: pathlib.Path) -> pathlib.Path:
        path = self._path / path
        if path.exists():
            raise RuntimeError(f'path already exists: "{path}"')

        path.mkdir()
        return path

    @staticmethod
    def rmtree(path: pathlib.Path | str):
        path = pathlib.Path(path)
        if path.exists():
            shutil.rmtree(path)

    def touch(self, path: pathlib.Path) -> pathlib.Path:
        if path in self._mappings:
            raise RuntimeError(f'file already exists: "{path}"')

        new_file = self._path / path
        new_file.touch()

        self._mappings[path] = new_file
        return new_file

    def _path_belongs_to_wd(self, path: pathlib.Path, check_if_exists: bool = True) -> bool:
        try:
            path.relative_to(self._path)
            return not check_if_exists or path.exists()
        except ValueError:
            return False

    def __getitem__(self, item: pathlib.Path | str) -> pathlib.Path:
        value = self.get(pathlib.Path(item))
        if value:
            return value

        raise KeyError(f'no such file "{item}"')

    def __contains__(self, item: pathlib.Path | str) -> bool:
        item = pathlib.Path(item)
        if self._path_belongs_to_wd(item):
            return True

        item, _ = self._parse_uri(item)
        return item in self._mappings

    def get(self, item: pathlib.Path | str, default=None):
        item = pathlib.Path(item)
        if self._path_belongs_to_wd(item):
            return item.resolve()

        item, grp = self._parse_uri(item)
        value = self._mappings.get(item, default)
        if value and grp:
            return pathlib.Path(f"{value}::{grp}")
        return value

    def get_staged_file_names(self) -> Dict[pathlib.Path, pathlib.Path]:
        return dict(sorted(self._mappings.items()))

    @property
    def name(self):
        return self._path

    def cleanup(self):
        if not self._delete or not self._path.exists():
            return

        def error_handler(_, path, excinfo):
            logging.warning(f'failed to delete "{path}": {excinfo}')

        shutil.rmtree(self._path, onexc=error_handler)


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


def _preprocess_plan(
    plan: Mapping[str, Any],
    wd: WorkingDirectory,
    fields_to_strip: List[str] | None = None,
) -> Tuple[bool, Dict[str, Any]]:
    skip = not _check_if_test_should_run(plan)
    digest = _hash_plan(plan, wd.name)
    plan = _strip_fields_from_config(plan, fields_to_strip)

    assert "id" not in plan
    plan["id"] = digest
    return skip, plan


def _get_uri(config: Dict[str, Any], fmt: str | None = None) -> pathlib.Path:
    if "files" not in config or len(config["files"]) == 0:
        raise ValueError("unable to fetch uri from config")

    if fmt is None:
        return pathlib.Path(config["files"][0]["uri"])

    for c in config["files"]:
        if c["format"] == fmt:
            return pathlib.Path(c["uri"])

    raise ValueError(f'unable to fetch uri with format "{fmt}" from config')


def _hash_plan(plan: Mapping[str, Any], tmpdir: pathlib.Path, algorithm: str = "sha256") -> str:
    tmpdir = str(tmpdir)

    def normalize_uri(uri: str) -> str:
        path, _, grp = uri.partition("::")
        path = pathlib.Path(path).name
        if grp:
            return f"{path}::{grp}"

        return path

    def strip_tmpdir(obj) -> Dict:
        for key, value in obj.items():
            if isinstance(value, dict) or isinstance(value, immutabledict):
                obj[key] = strip_tmpdir(value)
            elif isinstance(value, list):
                for i, x in enumerate(value):
                    if isinstance(x, str) and tmpdir in x:
                        value[i] = normalize_uri(x)
                obj[key] = value
            elif isinstance(value, str) and tmpdir in value:
                obj[key] = normalize_uri(value)

        return dict(sorted(obj.items()))

    serialized_plan = json.dumps(strip_tmpdir(dict(plan))).encode("utf-8")
    h = hashlib.new(algorithm)
    h.update(serialized_plan)

    return h.hexdigest()
