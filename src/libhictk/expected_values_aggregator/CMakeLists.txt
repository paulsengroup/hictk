# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

add_library(hictk_internal_expected_values_aggregator INTERFACE)
add_library(hictk::expected_values_aggregator ALIAS hictk_internal_expected_values_aggregator)

target_sources(
  hictk_internal_expected_values_aggregator
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_expected_values_aggregator
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_expected_values_aggregator
  INTERFACE
    hictk::common
    hictk::chromosome
    hictk::pixel
    hictk::reference
    hictk::transformers
)

target_link_system_libraries(
  hictk_internal_expected_values_aggregator
  INTERFACE
  phmap
  spdlog::spdlog_header_only
)
