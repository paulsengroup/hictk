# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(HDF5 QUIET COMPONENTS C)
if(NOT HDF5_FOUND)
  find_package(
    HDF5
    CONFIG
    QUIET
    REQUIRED
    COMPONENTS C)
endif()

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 QUIET)
  if(NOT Eigen3_FOUND)
    find_package(
      Eigen3
      CONFIG
      QUIET
      REQUIRED)
  endif()
endif()

add_library(cooler INTERFACE)
add_library(hictk::cooler ALIAS cooler)
target_sources(
  cooler
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
            ${CMAKE_CURRENT_SOURCE_DIR}/uri_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_equal_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/validation_impl.hpp)

target_include_directories(cooler INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(
  cooler
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
  cooler
  INTERFACE
  $<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>
  HDF5::HDF5
  std::filesystem
  _hictk_fast_float_tgt
  _hictk_fmt_tgt
  _hictk_highfive_tgt
  _hictk_phmap_tgt
  _hictk_spdlog_tgt)

install(TARGETS cooler INCLUDES)
