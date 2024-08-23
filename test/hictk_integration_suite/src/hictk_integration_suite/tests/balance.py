# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import is_hic


class _HictkBalanceCli(HictkTestHarness):
    def __init__(
        self,
        hictk_exec: pathlib.Path | str,
        algorithm: str,
        cwd: pathlib.Path | str | None = None,
        tmpdir: pathlib.Path | str | None = None,
    ):
        super().__init__(hictk_exec, cwd, tmpdir)

        if algorithm not in {"ice", "scale", "vc"}:
            raise NotImplementedError

        self._algorithm = algorithm

    def __repr__(self) -> str:
        return f"hictk-balance-{self._algorithm}-cli"

    def _validate(self, expect_failure: bool):
        if expect_failure:
            if len(self.stderr()) == 0:
                self._failures["missing help message"] = ""
            if len(self.stdout()) != 0:
                self._failures["unexpected output on stdout"] = self.stdout(500).strip()
            if self.returncode == 0:
                self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"
            return

        if len(self.stdout()) == 0:
            self._failures["missing help message"] = ""
        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        elif not expect_failure and self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"


class _HictkBalance(HictkTestHarness):
    def __init__(
        self,
        hictk_exec: pathlib.Path | str,
        algorithm: str,
        cwd: pathlib.Path | str | None = None,
        tmpdir: pathlib.Path | str | None = None,
    ):
        super().__init__(hictk_exec, cwd, tmpdir)

        if algorithm not in {"ice", "scale", "vc"}:
            raise NotImplementedError

        self._algorithm = algorithm

    def __repr__(self) -> str:
        return f"hictk-balance-{self._algorithm}"

    def _generate_normalization_name(self, mode: str):
        return f"{self._algorithm.upper()}_{mode.upper()}"

    def _validate(
        self,
        test_file: str,
        reference_file: str,
        resolutions: List[int],
        expect_failure: bool,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"

        mode = self._get_hictk_keyword_option("--mode", "gw")
        if mode not in {"cis", "trans", "gw"}:
            raise NotImplementedError

        normalization_name = self._generate_normalization_name(mode)

        if len(resolutions) == 0:
            resolutions = self._fetch_table(reference_file, table="resolutions")[1]
        for res in resolutions:
            _, expected = self._fetch_table(reference_file, resolution=res, table="normalizations")
            assert normalization_name in expected
            _, found = self._fetch_table(test_file, resolution=res, table="normalizations")

            if normalization_name not in found:
                self._failures[f"missing normalization {normalization_name}"] = ""
                continue

            # TODO remove check once hictkpy.File().weights() is available
            # https://github.com/paulsengroup/hictkpy/pull/49
            if not is_hic(reference_file) and not is_hic(test_file):
                _, expected = self._fetch_table(reference_file, resolution=res, table="weights")
                _, found = self._fetch_table(test_file, resolution=res, table="weights")
                self._failures |= validators.compare_weights(expected, found)

    def run(  # noqa
        self,
        args: List[str],
        file_format: str,
        variable_bin_size: bool,
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        expect_failure: bool = False,
        title: str | None = None,
        id: str | None = None,  # noqa
    ) -> Dict[str, Any]:
        if title is None:
            title = str(self)

        self.clear()
        self._id = id
        self._title = title
        self._args = args
        self._expect_failure = expect_failure

        t0 = timer()
        self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        t1 = timer()
        self._validate(
            file_format=file_format,
            variable_bin_size=variable_bin_size,
            expect_failure=expect_failure,
        )
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()


class HictkBalanceICE(_HictkBalance):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "ice"
        super().__init__(*args, **kwargs)


class HictkBalanceSCALE(_HictkBalance):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "scale"
        super().__init__(**kwargs)


class HictkBalanceVC(_HictkBalance):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "vc"
        super().__init__(**kwargs)


class HictkBalanceICECli(_HictkBalanceCli):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "ice"
        super().__init__(**kwargs)


class HictkBalanceSCALECli(_HictkBalanceCli):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "scale"
        super().__init__(**kwargs)


class HictkBalanceVCCli(_HictkBalanceCli):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "vc"
        super().__init__(**kwargs)
