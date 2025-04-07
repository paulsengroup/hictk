# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_pixel INTERFACE)
add_library(hictk::pixel ALIAS hictk_internal_pixel)

target_sources(
  hictk_internal_pixel
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/pixel.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/pixel_impl.hpp"
)

target_include_directories(
  hictk_internal_pixel
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_pixel
  INTERFACE
    hictk::bin_table
    hictk::chromosome
    hictk::format
    hictk::hash
    hictk::type_traits
)
