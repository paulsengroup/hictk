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

add_library(hictk_cooler INTERFACE)
add_library(hictk::cooler ALIAS hictk_cooler)
target_sources(
  hictk_cooler
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/attribute_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_accessors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_iterator_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_read_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_write_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_accessors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_read_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_standard_attr_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_validation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_write_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/index_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/multires_cooler_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/uri_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_copy_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_equal_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_merge_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/validation_impl.hpp)

target_include_directories(hictk_cooler INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(
  hictk_cooler
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::balancing
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::fmt
            hictk::genomic_interval
            hictk::pixel
            hictk::reference
            hictk::transformers
            hictk::variant)

target_link_system_libraries(
  hictk_cooler
  INTERFACE
  $<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>
  FastFloat::fast_float
  fmt::fmt
  HDF5::HDF5
  HighFive
  phmap
  spdlog::spdlog
  std::filesystem)

install(TARGETS hictk_cooler FILE_SET HEADERS EXCLUDE_FROM_ALL)
