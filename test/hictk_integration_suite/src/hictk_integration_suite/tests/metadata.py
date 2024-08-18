# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from hictk_integration_suite import validators
from hictk_integration_suite.runners.hictk import HictkTestHarness


class HictkMetadataCli(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-metadata-cli"

    def _validate(self, expect_failure: bool):
        if expect_failure:
            if len(self.stderr()) == 0:
                self.failures_["missing help message"] = ""
            if len(self.stdout()) != 0:
                self.failures_["unexpected output on stdout"] = self.stdout(100)
            if expect_failure and self.returncode() == 0:
                self.failures_["unexpected return code"] = f"expected non-zero, found {self.returncode()}"
            return

        if len(self.stdout()) == 0:
            self.failures_["missing help message"] = ""
        if len(self.stderr()) != 0:
            self.failures_["unexpected output on stderr"] = self.stderr(100)
        elif not expect_failure and self.returncode() != 0:
            self.failures_["unexpected return code"] = f"expected zero, found {self.returncode()}"


class HictkMetadata(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-metadata"

    def _validate(self, expect_failure: bool):
        if not expect_failure and len(self.stderr()) != 0:
            self.failures_["unexpected output on stderr"] = self.stderr(100)
        if expect_failure:
            raise NotImplemented()
        if not expect_failure and self.returncode() != 0:
            self.failures_["unexpected return code"] = f"expected zero, found {self.returncode()}"

        i = self.args().index("--output-format") + 1
        output_fmt = self.args()[i]
        if output_fmt not in {"json", "toml", "yaml"}:
            raise NotImplemented()

        payload = "".join(self.stdout())
        try:
            if output_fmt == "json":
                validators.metadata.json(payload)
            elif output_fmt == "toml":
                validators.metadata.toml(payload)
            elif output_fmt == "yaml":
                validators.metadata.yaml(payload)
        except (RuntimeError, ValueError) as e:
            self.failures_[f"output is not valid {output_fmt}"] = str(e)
