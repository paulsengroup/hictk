# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_common_tmpdir INTERFACE)
add_library(hictk::tmpdir ALIAS hictk_internal_common_tmpdir)

target_sources(
  hictk_internal_common_tmpdir
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/tmpdir.hpp"
)

target_include_directories(
  hictk_internal_common_tmpdir
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_common_tmpdir
  INTERFACE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
)
