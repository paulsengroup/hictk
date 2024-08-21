# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import inspect
import pathlib
import stat
import subprocess as sp
import sys
import textwrap

from hictk_integration_suite.runners.hictk import HictkTestHarness


class HictkTestHarnessImpl(HictkTestHarness):
    def _validate(self, expect_failure: bool):
        if expect_failure:
            if self._returncode == 0:
                self._failures[f"expected non-zero returncode, found {self._returncode}"] = ""
            if len(self._stdout) != 0:
                self._failures["unexpected output to stdout"] = ""
            if len(self._stderr) == 0:
                self._failures["missing output to stderr"] = ""
            return

        if self._returncode != 0:
            self._failures[f"expected zero returncode, found {self._returncode}"] = ""
        if len(self._stdout) == 0:
            self._failures["missing output to stdout"] = ""
        if len(self._stderr) != 0:
            self._failures["unexpected output to stderr"] = ""


class TestClass:
    @staticmethod
    def _mock_hictk(argv):
        import os
        import sys

        if len(argv) != 2:
            print("provide exactly one argument", file=sys.stderr)
            sys.exit(1)

        if "HICTK_FAIL" in os.environ:
            print("failing because HICTK_FAIL was found in the env variables", file=sys.stderr)
            sys.exit(1)

        exit_code = int(sys.argv[1])
        if exit_code == 0:
            print("All good!", file=sys.stdout)
            sys.exit(0)

        print("Something is wrong!", file=sys.stderr)
        sys.exit(exit_code)

    def _create_executable(self, dir: pathlib.Path, name: str = "hictk") -> pathlib.Path:  # noqa
        f = dir / name

        fx = textwrap.dedent(inspect.getsource(self._mock_hictk)).replace("@staticmethod", "")
        f.write_text(f"{fx}\n" "import sys\n" "_mock_hictk(sys.argv)\n")

        f.chmod(stat.S_IRUSR | stat.S_IXUSR)

        proc = sp.run([sys.executable, str(f), "0"], stdin=sp.DEVNULL, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        if proc.returncode != 0:
            raise RuntimeError("test script is broken")

        proc = sp.run([sys.executable, str(f)], stdin=sp.DEVNULL, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        if proc.returncode != 1:
            raise RuntimeError("test script is broken")

        proc = sp.run([sys.executable, str(f), "2"], stdin=sp.DEVNULL, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        if proc.returncode != 2:
            raise RuntimeError("test script is broken")

        return f

    def test_hictk_test_harness_ctor(self, tmp_path: pathlib.Path):
        exe = self._create_executable(tmp_path)
        (tmp_path / "tmp").mkdir()
        HictkTestHarnessImpl(exe, tmp_path, tmp_path / "tmp")

    def test_hictk_test_harness_run_success(self, tmp_path: pathlib.Path):
        exe = self._create_executable(tmp_path)
        (tmp_path / "tmp").mkdir()
        test = HictkTestHarnessImpl(sys.executable, tmp_path, tmp_path / "tmp")
        status = test.run([exe, "0"], timeout=1, env_variables={})
        assert test.ok()
        assert test.returncode == 0
        assert test.stderr() == ""
        assert test.stdout() == "All good!\n"
        assert test.stdout(3) == "All\n -- truncated"
        assert test.args == [sys.executable, str(exe), "0"]
        assert len(test.failures) == 0
        assert status["status"] == "PASS"

    def test_hictk_test_harness_run_failure(self, tmp_path: pathlib.Path):
        exe = self._create_executable(tmp_path)
        (tmp_path / "tmp").mkdir()
        test = HictkTestHarnessImpl(sys.executable, tmp_path, tmp_path / "tmp")
        status = test.run([exe, "10"], timeout=1, env_variables={})
        assert not test.ok()
        assert test.returncode == 10
        assert test.stderr() == "Something is wrong!\n"
        assert test.stderr(3) == "Som\n -- truncated"
        assert test.stdout() == ""

        assert test.args == [sys.executable, str(exe), "10"]
        assert len(test.failures) == 3
        assert status["status"] == "FAIL"

    def test_hictk_test_harness_run_failure_env(self, tmp_path: pathlib.Path):
        exe = self._create_executable(tmp_path)
        (tmp_path / "tmp").mkdir()
        test = HictkTestHarnessImpl(sys.executable, tmp_path, tmp_path / "tmp")
        status = test.run([exe, "0"], timeout=1, env_variables={"HICTK_FAIL": ""})
        assert status["status"] == "FAIL"

    def test_hictk_test_harness_run_success_w_expected_failure(self, tmp_path: pathlib.Path):
        exe = self._create_executable(tmp_path)
        (tmp_path / "tmp").mkdir()
        test = HictkTestHarnessImpl(sys.executable, tmp_path, tmp_path / "tmp")
        status = test.run([exe, "10"], timeout=1, env_variables={}, expect_failure=True)
        assert not test.ok()
        assert test.returncode == 10
        assert len(test.failures) == 0
        assert status["status"] == "PASS"
