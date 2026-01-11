# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_file STATIC)
add_library(hictk::file ALIAS hictk_internal_file)

target_sources(
  hictk_internal_file
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/file.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/multires_file.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/file_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/multires_file_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/utils_merge_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/multires_file.cpp"
)

target_include_directories(
  hictk_internal_file
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_file
  PUBLIC
    hictk::balancing
    hictk::bin_table
    hictk::cooler
    hictk::genomic_interval
    hictk::hic
    hictk::pixel
    hictk::reference
    hictk::tmpdir
)

target_link_system_libraries(hictk_internal_file PUBLIC fmt::fmt-header-only)
