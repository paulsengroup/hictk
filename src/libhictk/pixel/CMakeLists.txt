# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_pixel STATIC)
add_library(hictk::pixel ALIAS hictk_internal_pixel)

target_sources(
  hictk_internal_pixel
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/pixel.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/pixel_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel.cpp"
)

target_include_directories(
  hictk_internal_pixel
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_pixel
  PUBLIC
    hictk::bin_table
    hictk::chromosome
    hictk::format
    hictk::hash
    hictk::type_traits
    hictk::numeric
  PRIVATE
    hictk::genomic_interval
)

target_link_system_libraries(hictk_internal_pixel PRIVATE fmt::fmt-header-only)
