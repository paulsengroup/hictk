# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(HDF5 REQUIRED QUIET COMPONENTS C)
find_package(HighFive REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_cooler STATIC)
add_library(hictk::cooler ALIAS hictk_internal_cooler)
target_sources(
  hictk_internal_cooler
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
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
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_merge_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/utils_rename_chroms_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/cooler/impl/validation_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_accessors.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_read.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_write.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_accessors.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_read.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_standard_attr.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_write.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/group.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/index.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/multires_cooler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/singlecell_cooler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/uri.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_copy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_equal.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/validation.cpp"
)

target_include_directories(
  hictk_internal_cooler
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_cooler
  PUBLIC
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
  PUBLIC
    fmt::fmt-header-only
    HDF5::HDF5
    "HighFive$<$<VERSION_GREATER_EQUAL:${HighFive_VERSION},3>:::HighFive>"
    phmap
    spdlog::spdlog_header_only
)
