# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(HighFive REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_cooler INTERFACE)
add_library(hictk::cooler ALIAS hictk_internal_cooler)
target_sources(
  hictk_internal_cooler
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_cooler
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_cooler
  INTERFACE
    hictk::balancing
    hictk::bin_table
    hictk::chromosome
    hictk::common
    hictk::format
    hictk::genomic_interval
    hictk::pixel
    hictk::reference
    hictk::transformers
    hictk::type_pretty_printer
    hictk::type_traits
    hictk::variant
)

target_link_system_libraries(
  hictk_internal_cooler
  INTERFACE
  FastFloat::fast_float
  fmt::fmt-header-only
  HDF5::HDF5
  HighFive
  phmap
  spdlog::spdlog_header_only
)
