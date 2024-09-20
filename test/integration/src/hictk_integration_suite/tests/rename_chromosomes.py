# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import logging
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from hictk_integration_suite.runners.hictk import HictkTestHarness
from hictk_integration_suite.validators.file_formats import (
    is_cooler,
    is_multires,
    is_scool,
)

from .cli import HictkCli


class HictkRenameChromosomesCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-rename-chromosomes-cli"


class HictkRenameChromosomes(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-rename-chromosomes"

    def _validate_cool(self, uri: pathlib.Path | str, action: str):
        _, found = self._fetch_table(uri, table="chroms")
        assert found is not None

        if action == "add_prefix":
            for chrom in found:
                num_prefix = chrom.count("chr")
                if num_prefix != 1 and num_prefix != 2:
                    self._failures[f'unexpected prefix for "{chrom}"'] = (
                        f'expected one or two intances of "chr" prefix, found {num_prefix}'
                    )
        elif action == "remove_prefix":
            for chrom in found:
                num_prefix = chrom.count("chr")
                if num_prefix != 0:
                    self._failures[f'unexpected prefix for "{chrom}"'] = (
                        f'expected zero intance of "chr" prefix, found {num_prefix}'
                    )
        else:
            raise NotImplementedError

    def _validate_mcool(self, file: pathlib.Path | str, action: str):
        resolutions = self._fetch_table(file, table="resolutions")
        assert resolutions is not None

        for res in resolutions:
            self._validate_cool(f"{file}::/resolutions/{res}", action)

    def _validate_scool(self, file: pathlib.Path | str, action: str):
        cells = self._fetch_table(file, table="cells")
        assert cells is not None and len(cells) != 0

        return self._validate_cool(f"{file}::/cells/{cells[0]}", action)

    def _validate_cool_name_mappings(self, uri: pathlib.Path | str, name_mappings: Dict[str, str]):
        _, found = self._fetch_table(uri, table="chroms")
        assert found is not None

        for src, dest in name_mappings.items():
            if src in found:
                self._failures[f'found chromosome with original name "{src}"'] = f'should\'ve been renamed to "{dest}"'
            elif dest not in found:
                self._failures[f'missing chromosome "{dest}"'] = ""

    def _validate_mcool_name_mappings(self, file: pathlib.Path | str, name_mappings: Dict[str, str]):
        resolutions = self._fetch_table(file, table="resolutions")
        assert resolutions is not None

        for res in resolutions:
            self._validate_cool_name_mappings(f"{file}::/resolutions/{res}", name_mappings)

    def _validate_scool_name_mappings(self, file: pathlib.Path | str, name_mappings: Dict[str, str]):
        cells = self._fetch_table(file, table="cells")
        assert cells is not None and len(cells) != 0

        return self._validate_cool_name_mappings(f"{file}::/cells/{cells[0]}", name_mappings)

    def _validate(self, test_file: str, expect_failure: bool, name_mappings: str | None):
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        action = None
        if self._get_hictk_flag_value("--add-chr-prefix"):
            action = "add_prefix"
        elif self._get_hictk_flag_value("--remove-chr-prefix"):
            action = "remove_prefix"

        if name_mappings is None:
            if is_cooler(test_file):
                self._validate_cool(test_file, action)
            elif is_multires(test_file):
                self._validate_mcool(test_file, action)
            elif is_scool(test_file):
                self._validate_scool(test_file, action)
            else:
                raise NotImplementedError
            return

        assert pathlib.Path(name_mappings).is_file()
        name_mappings = pd.read_table(name_mappings, names=["src", "dest"]).set_index("src")["dest"].to_dict()

        if is_cooler(test_file):
            self._validate_cool_name_mappings(test_file, name_mappings)
        elif is_multires(test_file):
            self._validate_mcool_name_mappings(test_file, name_mappings)
        elif is_scool(test_file):
            self._validate_scool_name_mappings(test_file, name_mappings)
        else:
            raise NotImplementedError

    def run(  # noqa
        self,
        args: List[str],
        test_file: str,
        name_mappings: str | None = None,
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
        try:
            self._run_hictk(args, timeout=timeout, env_variables=env_variables)
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(
                test_file=test_file,
                name_mappings=name_mappings,
                expect_failure=expect_failure,
            )
        except:  # noqa
            logging.error(f"failed to validate output produced by {args}")
            raise
        t2 = timer()

        self._hictk_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
