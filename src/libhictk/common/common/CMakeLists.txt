# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_common INTERFACE)
add_library(hictk::common ALIAS hictk_internal_common)

target_sources(
  hictk_internal_common
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/common.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/suppress_warnings.hpp"
)

target_include_directories(
  hictk_internal_common
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
