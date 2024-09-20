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
import tempfile
from typing import Any, Dict, List, Mapping, Tuple

from immutabledict import immutabledict

from hictk_integration_suite.common import URI


def _make_file_read_only(path: pathlib.Path | str):
    mode = stat.S_IRUSR
    if os.access(path, os.X_OK):
        mode |= stat.S_IXUSR
    pathlib.Path(path).chmod(mode)


def _make_file_writeable(path: pathlib.Path | str):
    mode = stat.S_IRUSR | stat.S_IWUSR
    if os.access(path, os.X_OK):
        mode |= stat.S_IXUSR
    pathlib.Path(path).chmod(mode)


def _argument_map_to_list(args_map: Dict[str, Any]) -> List[str]:
    args = []
    for k, v in args_map.items():
        k = "--" + (str(k).removeprefix("--"))
        if v is None:
            args.append(k)
        elif isinstance(v, str) and len(v) == 0:
            args.append(k)
        elif isinstance(v, list):
            args.append(k)
            args.extend((str(x) for x in v))
        else:
            args.extend((k, str(v)))

    return args


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
    def _make_read_only(path: URI):
        _make_file_read_only(path.path)

    @staticmethod
    def _make_writeable(path: URI):
        _make_file_writeable(path.path)

    def stage_file(
        self,
        src: URI | pathlib.Path | str,
        make_read_only: bool = True,
        exists_ok: bool = False,
    ) -> pathlib.Path:

        src = URI(src, False)
        if not src.path.exists():
            raise RuntimeError(f'source file "{src.path}" does not exist')

        dest_dir = self._path / "staged_files"
        if src.group is None:
            dest = URI(dest_dir / src.path.name)
        else:
            dest = URI(dest_dir / f"{src.path.name}::{src.group}")

        if dest.path.exists():
            if not exists_ok:
                raise RuntimeError(f'refusing to overwrite file "{dest.path}"')

            logging.debug(f'file "{src.path}" was already staged')
        else:
            logging.debug(f'staging file "{src.path}"...')
            dest_dir.mkdir(exist_ok=True)
            shutil.copy2(src.path, dest.path)
            if make_read_only:
                self._make_read_only(dest)

        self._mappings[src] = dest
        self._mappings[dest] = dest
        if src.group is not None:
            assert dest.group is not None
            dest_file = URI(dest.path)
            self._mappings[URI(src.path)] = dest_file
            self._mappings[dest_file] = dest_file

        logging.debug(f'URI "{src}" successfully staged (dest="{dest}")')
        return dest

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

    def __getitem__(self, item: URI | pathlib.Path | str) -> URI:
        value = self.get(URI(item))
        if value is not None:
            return value

        raise KeyError(f'no such file "{item}"')

    def __contains__(self, item: URI | pathlib.Path | str) -> bool:
        item = URI(item)
        if self._path_belongs_to_wd(item.path) and item.path.exists():
            return True

        return item in self._mappings

    def get(self, item: URI | pathlib.Path | str, default=None) -> URI | None:
        item = URI(item)
        if self._path_belongs_to_wd(item.path) and item.path.exists():
            return item

        return self._mappings.get(item, default)

    def get_staged_file_names(self) -> Dict[str, str]:
        return {str(k): str(v) for k, v in sorted(self._mappings.items())}

    @property
    def name(self):
        return self._path

    def cleanup(self):
        if not self._delete or not self._path.exists():
            return

        # Make files writeable by the current user
        for dirpath, dirnames, filenames in self._path.walk(follow_symlinks=False):
            for name in dirnames:
                try:
                    path = dirpath / name
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                except:  # noqa
                    pass
            for name in filenames:
                try:
                    path = dirpath / name
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                except:  # noqa
                    pass

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

    def strip_tmpdir(obj) -> Dict:
        for key, value in obj.items():
            if isinstance(value, dict) or isinstance(value, immutabledict):
                obj[key] = strip_tmpdir(value)
            elif isinstance(value, list):
                for i, x in enumerate(value):
                    if isinstance(x, str) and tmpdir in x:
                        value[i] = str(URI(x))
                obj[key] = value
            elif isinstance(value, str) and tmpdir in value:
                obj[key] = str(URI(value))

        return dict(sorted(obj.items()))

    serialized_plan = json.dumps(strip_tmpdir(dict(plan))).encode("utf-8")
    h = hashlib.new(algorithm)
    h.update(serialized_plan)

    return h.hexdigest()
