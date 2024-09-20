# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
import stat
import tempfile
from typing import Tuple

import pytest

from hictk_integration_suite.cli.common import URI, WorkingDirectory


class TestClass:
    @staticmethod
    def _mkdtemp():
        return tempfile.TemporaryDirectory(prefix="hictk-integration-suite-")

    @staticmethod
    def _create_test_files(tmpdir: str) -> Tuple[pathlib.Path, pathlib.Path]:
        plain_file = pathlib.Path(tmpdir) / "plain.txt"
        exec_file = pathlib.Path(tmpdir) / "script.sh"

        plain_file.write_text("foo\n")
        plain_file.chmod(stat.S_IRUSR | stat.S_IWUSR)

        exec_file.write_text("#!/usr/bin/env bash\necho Hi!")
        exec_file.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

        return plain_file, exec_file

    @staticmethod
    def _is_executable(f: URI) -> bool:
        return os.access(f.path, os.X_OK)

    @staticmethod
    def _attempt_write(f: URI) -> bool:
        assert f.path.is_file()
        f.path.write_text("foo")
        return True

    def test_default_ctor(self):
        wd = WorkingDirectory()
        path = wd.name
        assert path.is_dir()
        wd.cleanup()
        assert not path.exists()

        if path.exists():
            path.rmdir()

    def test_existing_custom_tmpdir_ctor(self):
        with self._mkdtemp() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)

            with pytest.raises(RuntimeError):
                WorkingDirectory(tmpdir)

            assert tmpdir.is_dir()
            WorkingDirectory(tmpdir / "foo")

    def test_non_existing_custom_tmpdir_ctor(self):
        with self._mkdtemp() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            new_tmpdir = tmpdir / "foo"
            with WorkingDirectory(new_tmpdir, delete=True) as wd:
                assert wd.name == new_tmpdir
                assert new_tmpdir.is_dir()
            assert not new_tmpdir.exists()

            new_tmpdir = tmpdir / "bar"
            with WorkingDirectory(new_tmpdir, delete=False) as wd:
                assert wd.name == new_tmpdir
                assert new_tmpdir.is_dir()
            assert new_tmpdir.is_dir()

    def test_file_staging_plain_file_ro(self):
        with self._mkdtemp() as tmpdir:
            plain_file, _ = self._create_test_files(tmpdir)
            with WorkingDirectory() as wd:
                assert plain_file not in wd
                wd.stage_file(plain_file)
                f = wd.get(plain_file)
                assert f is not None
                if f:
                    assert not self._is_executable(f)
                    with pytest.raises(Exception):
                        self._attempt_write(f)

    def test_file_staging_exec_file_ro(self):
        with self._mkdtemp() as tmpdir:
            _, exec_file = self._create_test_files(tmpdir)
            with WorkingDirectory() as wd:
                assert exec_file not in wd
                wd.stage_file(exec_file)
                f = wd.get(exec_file)
                assert f is not None
                if f:
                    assert self._is_executable(f)
                    with pytest.raises(Exception):
                        self._attempt_write(f)

    def test_file_staging_plain_file_rw(self):
        with self._mkdtemp() as tmpdir:
            plain_file, _ = self._create_test_files(tmpdir)
            with WorkingDirectory() as wd:
                assert plain_file not in wd
                wd.stage_file(plain_file, make_read_only=False)
                f = wd.get(plain_file)
                assert f is not None
                if f:
                    assert not self._is_executable(f)
                    assert self._attempt_write(f)

    def test_file_staging_exec_file_rw(self):
        with self._mkdtemp() as tmpdir:
            _, exec_file = self._create_test_files(tmpdir)
            with WorkingDirectory() as wd:
                assert exec_file not in wd
                wd.stage_file(exec_file, make_read_only=False)
                f = wd.get(exec_file)
                assert f is not None
                if f:
                    assert self._is_executable(f)
                    assert self._attempt_write(f)

    def test_file_staging_uri(self):
        with self._mkdtemp() as tmpdir:
            plain_file, _ = self._create_test_files(tmpdir)
            uri = f"{plain_file}::/group/foobar"
            with WorkingDirectory() as wd:
                assert plain_file not in wd
                assert uri not in wd
                wd.stage_file(uri)
                assert wd.get(plain_file) is not None
                assert wd.get(uri) is not None

    def test_mkdtemp(self):
        with WorkingDirectory() as wd:
            assert wd.mkdtemp().is_dir()

    def test_mkdir(self):
        with WorkingDirectory() as wd:
            assert wd.mkdir("foo").is_dir()
            with pytest.raises(RuntimeError):
                wd.mkdir("foo")

    def test_touch(self):
        with WorkingDirectory() as wd:
            assert wd.touch("foo").is_file()
            wd.mkdir("dir")
            assert wd.touch("dir/foo").is_file()
            with pytest.raises(RuntimeError):
                wd.touch("foo")

    def test_get_staged_file_names(self):
        with WorkingDirectory() as wd:
            wd.touch("foo").is_file()
            wd.mkdir("dir")

            files = wd.get_staged_file_names()
            assert "foo" in files
