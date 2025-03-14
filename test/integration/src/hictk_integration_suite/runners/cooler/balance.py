# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import warnings
from typing import List, Tuple

import cooler
import numpy.typing as npt
import structlog
from hictk_integration_suite.validators.file_formats import is_cooler, is_multires


class CoolerBalance:
    def __init__(
        self,
        uri: pathlib.Path,
        resolution: int | None,
    ):
        self._clr = self._open(str(uri), resolution)

    @staticmethod
    def _open(uri: str, resolution: int | None) -> cooler.Cooler:
        if is_cooler(uri):
            clr = cooler.Cooler(uri)
            if resolution is not None:
                assert clr.binsize == resolution
            return clr
        if is_multires(uri) and resolution is not None:
            return cooler.Cooler(f"{uri}::/resolutions/{resolution}")

        raise RuntimeError(f'URI "{uri}" does not point to a .[m]cool file')

    @staticmethod
    def _generate_normalization_name(mode: str) -> str:
        if mode == "cis":
            return "ICE"
        if mode == "trans":
            mode = "inter"
        return f"{mode.upper()}_ICE"

    def balance(self, mode: str, **kwargs) -> Tuple[str, npt.NDArray]:
        if mode == "cis":
            kwargs["cis_only"] = True
            kwargs["trans_only"] = False
        elif mode == "gw":
            kwargs["cis_only"] = False
            kwargs["trans_only"] = False
        elif mode == "trans":
            kwargs["cis_only"] = False
            kwargs["trans_only"] = True
        else:
            raise NotImplementedError

        name = self._generate_normalization_name(mode)

        if name in self._clr.bins():
            return name, self._clr.bins()[name][:].tonumpy()

        structlog.get_logger().debug(f"balancing cooler at URI {self._clr.uri} with {name}...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bias, _ = cooler.balance_cooler(self._clr, **kwargs)

        return name, bias
