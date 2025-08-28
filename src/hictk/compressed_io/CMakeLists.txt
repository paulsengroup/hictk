# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(LibArchive REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_compressed_io STATIC)

target_sources(hictk_internal_tools_compressed_io PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/compressed_reader.cpp")

target_include_directories(hictk_internal_tools_compressed_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_compressed_io
  PRIVATE
    hictk::project_options
    hictk::project_warnings
)

target_link_system_libraries(
  hictk_internal_tools_compressed_io
  PRIVATE
    fmt::fmt-header-only
    std::filesystem
  PUBLIC
    LibArchive::LibArchive
)
