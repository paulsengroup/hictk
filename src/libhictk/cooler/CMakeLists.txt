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
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/attribute.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/common.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/cooler.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/dataset.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/group.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/index.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/multires_cooler.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/pixel_selector.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/singlecell_cooler.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/uri.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/utils.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/validation.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/attribute_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/dataset_accessors_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/dataset_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/dataset_iterator_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/dataset_read_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/dataset_write_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_accessors_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_read_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_standard_attr_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_validation_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/file_write_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/index_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/multires_cooler_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/pixel_selector_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/singlecell_cooler_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/uri_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_copy_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_equal_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_merge_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_rename_chroms_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/validation_impl.hpp"
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
    hictk::string
    hictk::transformers
    hictk::type_pretty_printer
    hictk::type_traits
    hictk::variant
    hictk::version
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
