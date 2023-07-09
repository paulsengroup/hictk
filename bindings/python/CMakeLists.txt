# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FetchContent_Declare(
  _hictk_pybind11
  URL ${PROJECT_SOURCE_DIR}/external/pybind11-v2.10.3.tar.xz
  URL_HASH SHA256=1a0ab8cdfdab3a1e959f67a1f4b231e4dc483b545001051d3411981827031920
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_MakeAvailable(_hictk_pybind11)

add_subdirectory(src)
