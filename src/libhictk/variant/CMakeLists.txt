# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_variant INTERFACE)
add_library(hictk::variant ALIAS hictk_internal_variant)

target_sources(
  hictk_internal_variant
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/generic_variant.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/numeric_variant.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/variant_buff.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/variant_buff_impl.hpp"
)

target_include_directories(
  hictk_internal_variant
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
