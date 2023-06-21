# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(HighFive CONFIG QUIET)
if (NOT HighFive_FOUND)
  find_package(HighFive QUIET)
endif()

find_package(HDF5 CONFIG QUIET)
if (NOT HDF5_FOUND)
  find_package(HDF5 QUIET REQUIRED)
endif()

add_library(hictk INTERFACE)
add_library(hictk::cooler ALIAS hictk)
target_sources(
  hictk
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/attribute_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/balancing_impl.hpp
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

target_include_directories(hictk INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

# cmake-format: off
target_link_libraries(
  hictk
  INTERFACE hictk_project_warnings hictk_project_options hictk::bin_table hictk::chromosome hictk::common hictk::fmt hictk::genomic_interval hictk::pixel hictk::reference hictk::variant
)
# cmake-format: on

target_link_system_libraries(
  hictk
  INTERFACE
  FastFloat::fast_float
  fmt::fmt
  HDF5::HDF5
  HighFive
  std::filesystem)

install(TARGETS hictk INCLUDES)
