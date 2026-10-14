# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_genomic_interval STATIC)
add_library(hictk::genomic_interval ALIAS hictk_internal_genomic_interval)

target_sources(
  hictk_internal_genomic_interval
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/genomic_interval.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/genomic_interval_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/genomic_interval.cpp"
)

target_include_directories(
  hictk_internal_genomic_interval
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_genomic_interval
  PUBLIC
    hictk::chromosome
    hictk::reference
  PRIVATE
    hictk::hash
    hictk::numeric
)

target_link_system_libraries(hictk_internal_genomic_interval PRIVATE fmt::fmt-header-only)
