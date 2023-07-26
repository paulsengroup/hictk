# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(FastFloat CONFIG REQUIRED)
find_package(FMT CONFIG REQUIRED)
find_package(
  HDF5
  CONFIG
  REQUIRED
  QUIET
  COMPONENTS C)
find_package(HighFive CONFIG REQUIRED)
find_package(phmap CONFIG REQUIRED)
find_package(spdlog CONFIG REQUIRED)

if(HICTK_WITH_EIGEN)
  find_package(
    Eigen3
    CONFIG
    REQUIRED
    QUIET)
endif()

add_library(cooler INTERFACE)
add_library(hictk::cooler ALIAS cooler)
target_sources(
  cooler
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(cooler INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                            $<INSTALL_INTERFACE:include>)

target_link_libraries(
  cooler
  INTERFACE hictk::balancing
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::format
            hictk::genomic_interval
            hictk::pixel
            hictk::reference
            hictk::transformers
            hictk::variant)

target_link_system_libraries(
  cooler
  INTERFACE
  $<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>
  FastFloat::fast_float
  fmt::fmt
  HDF5::HDF5
  HighFive
  phmap
  spdlog::spdlog
  std::filesystem)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
