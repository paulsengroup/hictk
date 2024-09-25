# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import json
import logging
import tomllib
from timeit import default_timer as timer
from typing import Any, Dict, List

import yaml

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness

from .cli import HictkCli


class HictkMetadataCli(HictkCli):
    def __repr__(self) -> str:
        return "hictk-metadata-cli"


class HictkMetadata(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-metadata"

    def _validate_cool(self, data: Dict[str, Any], variable_bin_size: bool):
        format = data.get("format")  # noqa
        if format is None:
            self._failures["format attribute is missing"] = ""
        elif format != "HDF5::Cooler":
            self._failures["invalid format attribute"] = f'expected HDF5::Cooler, found "{format}"'

        bin_type = data.get("bin-type")
        if bin_type is None:
            self._failures["bin-type attribute is missing"] = ""
        elif variable_bin_size and bin_type != "variable":
            self._failures["invalid bin-type attribute"] = f'expected "variable", found "{bin_type}"'
        elif not variable_bin_size and bin_type != "fixed":
            self._failures["invalid bin-type attribute"] = f'expected "fixed", found "{bin_type}"'

    def _validate_mcool(self, data: Dict[str, Any]):
        format = data.get("format")  # noqa
        if format is None:
            self._failures["format attribute is missing"] = ""
        elif format != "HDF5::MCOOL":
            self._failures["invalid format attribute"] = f'expected HDF5::MCOOL, found "{format}"'

    def _validate_scool(self, data: Dict[str, Any]):
        format = data.get("format")  # noqa
        if format is None:
            self._failures["format attribute is missing"] = ""
        elif format != "HDF5::SCOOL":
            self._failures["invalid format attribute"] = f'expected HDF5::SCOOL, found "{format}"'

    def _validate_hic(self, data: Dict[str, Any]):
        format = data.get("format")  # noqa
        if format is None:
            self._failures["format attribute is missing"] = ""
        elif format != "HIC":
            self._failures["invalid format attribute"] = f'expected HIC, found "{format}"'

    def _validate(self, expect_failure: bool, file_format: str, variable_bin_size: bool):
        if expect_failure:
            raise NotImplementedError

        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"

        output_fmt = self._get_hictk_keyword_option("--output-format", "json")
        if output_fmt not in {"json", "toml", "yaml"}:
            raise NotImplementedError

        payload = "".join(self.stdout())
        data = None
        try:
            if output_fmt == "json":
                validators.metadata.json(payload)
                data = json.loads(payload)
            elif output_fmt == "toml":
                validators.metadata.toml(payload)
                data = tomllib.loads(payload)
            elif output_fmt == "yaml":
                validators.metadata.yaml(payload)
                data = yaml.safe_load(payload)
        except (RuntimeError, ValueError) as e:
            self._failures[f"output is not valid {output_fmt}"] = str(e)

        if data is None:
            return

        if file_format == "cool":
            self._validate_cool(data, variable_bin_size)
        elif file_format == "mcool":
            self._validate_mcool(data)
        elif file_format == "scool":
            self._validate_scool(data)
        elif file_format == "hic":
            self._validate_hic(data)
        else:
            raise NotImplementedError

    def run(  # noqa
        self,
        args: List[str],
        file_format: str,
        variable_bin_size: bool,
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        max_attempts: int = 1,
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
            self._run_hictk(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)
        except:  # noqa
            logging.error(f"failed to execute {args}")
            raise
        t1 = timer()
        try:
            self._validate(
                file_format=file_format,
                variable_bin_size=variable_bin_size,
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
