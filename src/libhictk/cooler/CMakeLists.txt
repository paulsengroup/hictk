# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat REQUIRED)
find_package(FMT REQUIRED)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(HighFive REQUIRED)
find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 REQUIRED QUIET)
endif()

add_library(cooler INTERFACE)
add_library(hictk::cooler ALIAS cooler)
target_sources(
  cooler
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  cooler
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  cooler
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
    hictk::variant
)

target_link_system_libraries(
  cooler
  INTERFACE
  "$<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>"
  FastFloat::fast_float
  fmt::fmt-header-only
  HDF5::HDF5
  HighFive
  phmap
  spdlog::spdlog_header_only
)
