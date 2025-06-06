# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_genomic_interval INTERFACE)
add_library(hictk::genomic_interval ALIAS hictk_internal_genomic_interval)

target_sources(
  hictk_internal_genomic_interval
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/genomic_interval.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/genomic_interval_impl.hpp"
)

target_include_directories(
  hictk_internal_genomic_interval
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_genomic_interval
  INTERFACE
    hictk::chromosome
    hictk::hash
    hictk::numeric
)
