# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(default_delete_libdeflate)
add_subdirectory(default_delete_zstd)
add_subdirectory(genomic_units)
add_subdirectory(type_pretty_printer)
add_subdirectory(type_traits)

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_common INTERFACE)
add_library(hictk::common ALIAS hictk_internal_common)
add_dependencies(hictk_internal_common _hictk_check_git)

target_sources(
  hictk_internal_common
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_common
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_common
  INTERFACE
    hictk::numeric
    fmt::fmt-header-only
    spdlog::spdlog_header_only
)
