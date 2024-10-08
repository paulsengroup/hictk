# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import warnings
from typing import List, Tuple

import cooler
import numpy.typing as npt
from hictk_integration_suite.validators.file_formats import is_cooler, is_multires


class CoolerCoarsen:
    def __init__(
        self,
        input_uri: pathlib.Path,
        base_resolution: int | None,
        output_uri: pathlib.Path,
        target_resolution: int,
    ):
        assert target_resolution > 0
        clr = self._open(str(input_uri), base_resolution)
        self._input_uri = clr.uri
        self._output_uri = output_uri

        if clr.binsize >= target_resolution or target_resolution % clr.binsize != 0:
            raise RuntimeError(
                f"CoolerCoarsen(): target_resolution should be a multiple of base resolution: found {target_resolution} and {base_resolution}, respectively"
            )

        self._factor = target_resolution // clr.binsize

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

    def coarsen(self, **kwargs) -> pathlib.Path:
        kwargs["factor"] = self._factor

        if "chunksize" not in kwargs:
            kwargs["chunksize"] = 10_000_000

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cooler.coarsen_cooler(self._input_uri, str(self._output_uri), **kwargs)

        return self._output_uri
