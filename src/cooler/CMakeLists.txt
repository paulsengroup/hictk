# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat CONFIG REQUIRED)
find_package(Filesystem REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(HDF5 CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)
find_package(tsl-hopscotch-map CONFIG REQUIRED)

add_library(coolerpp INTERFACE)
add_library(Coolerpp::Coolerpp ALIAS coolerpp)
target_sources(
  coolerpp
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/attribute_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/balancing_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/bin_table_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/chromosome_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_accessors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_read_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_standard_attr_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_validation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/coolerpp_write_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_accessors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_iterator_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_read_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/dataset_write_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/index_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/pixel_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/uri_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_equal_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/validation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/variant_buff_impl.hpp)

target_include_directories(coolerpp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

# cmake-format: off
target_link_libraries(
  coolerpp
  INTERFACE coolerpp_project_warnings coolerpp_project_options
)
# cmake-format: on

target_link_system_libraries(
  coolerpp
  INTERFACE
  FastFloat::fast_float
  fmt::fmt
  HDF5::HDF5
  HighFive
  tsl::hopscotch_map
  std::filesystem)

install(TARGETS coolerpp INCLUDES)
