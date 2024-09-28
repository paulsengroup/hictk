# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import os
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from hictk_integration_suite import validators
from hictk_integration_suite.runners.cooler import CoolerBalance
from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import is_hic

from .cli import HictkCli


class _HictkBalanceCli(HictkCli):
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

    def _generate_normalization_name(self, mode: str) -> str:
        if mode == "cis":
            return self._algorithm.upper()
        if mode == "trans":
            mode = "inter"
        return f"{mode.upper()}_{self._algorithm.upper()}"

    def _balance_clr(self, uri: str, resolution: int | None, mode: str) -> pd.DataFrame:
        name = self._generate_normalization_name(mode)
        return pd.DataFrame({name: CoolerBalance(uri, resolution).balance(mode, max_iters=500)[1]})

    def _validate(
        self,
        test_file: str,
        reference_file: str,
        expect_failure: bool,
        no_validate_weights: bool,
        resolutions: List[int] = None,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        mode = self._get_hictk_keyword_option("--mode", "gw")
        if mode not in {"cis", "trans", "gw"}:
            raise NotImplementedError

        normalization_name = self._generate_normalization_name(mode)

        if resolutions is None or len(resolutions) == 0:
            resolutions = self._fetch_table(reference_file, table="resolutions")[1]
        for res in resolutions:
            _, found = self._fetch_table(test_file, resolution=res, table="normalizations")

            if normalization_name not in found:
                self._failures[f"missing normalization {normalization_name}"] = ""
                continue

            if no_validate_weights:
                continue

            # TODO remove check once hictkpy.File().weights() is available
            # https://github.com/paulsengroup/hictkpy/pull/49
            if not is_hic(reference_file) and not is_hic(test_file):
                _, expected = self._fetch_table(reference_file, resolution=res, table="weights")
                if expected is None or normalization_name not in expected:
                    expected = self._balance_clr(reference_file, res, mode)

                _, found = self._fetch_table(test_file, resolution=res, table="weights")
                self._failures |= validators.compare_weights(
                    expected[[normalization_name]],
                    found[[normalization_name]],
                    atol=1.0,
                )

    def run(  # noqa
        self,
        args: List[str],
        reference_file: pathlib.Path | str,
        test_file: pathlib.Path | str,
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        max_attempts: int = 1,
        no_validate_weights: bool = False,
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

        if env_variables is None:
            env_variables = os.environ.copy()
        else:
            env_variables = dict(env_variables.copy())

        if "LLVM_PROFILE_FILE" in env_variables:
            env_variables["LLVM_PROFILE_FILE"] = env_variables["LLVM_PROFILE_FILE"].replace("%id", str(id))

        t0 = timer()
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)

        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(
                reference_file=reference_file,
                test_file=test_file,
                expect_failure=expect_failure,
                no_validate_weights=no_validate_weights,
            )
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()


class HictkBalanceICE(_HictkBalance):
    def __init__(self, **kwargs):
        kwargs["algorithm"] = "ice"
        super().__init__(**kwargs)


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
