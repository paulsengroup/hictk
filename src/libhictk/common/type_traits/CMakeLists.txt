# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_common_type_traits INTERFACE)
add_library(hictk::type_traits ALIAS hictk_internal_common_type_traits)

target_sources(
  hictk_internal_common_type_traits
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/type_traits.hpp"
)

target_include_directories(
  hictk_internal_common_type_traits
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
