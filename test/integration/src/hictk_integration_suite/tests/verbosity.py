# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from hictk_integration_suite.runners.hictk import HictkTestHarness


class HictkVerbosity(HictkTestHarness):
    def __repr__(self) -> str:
        return "hictk-verbosity"

    def _validate(
        self,
        *args,
        **kwargs,
    ):  # noqa
        if self.returncode != 0:
            self._failures["unexpected return code"] = (
                f"expected zero, found {self.returncode}; excerpt from stderr: {self.stderr(500).strip()}"
            )
            return
