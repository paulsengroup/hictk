# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(readerwriterqueue REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

add_library(hictk_internal_tools_convert STATIC)

target_sources(
  hictk_internal_tools_convert
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/convert.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cool_to_hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_to_cool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/common.hpp"
)

target_include_directories(hictk_internal_tools_convert PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_convert
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk::balancing
    hictk::bin_table
    hictk::chromosome
    hictk::cooler
    hictk::hic
    hictk::pixel
    hictk::reference
    hictk::string
    hictk::tmpdir
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_convert
  PRIVATE
    fmt::fmt-header-only
    readerwriterqueue::readerwriterqueue
    spdlog::spdlog_header_only
    std::filesystem
  PUBLIC
    Threads::Threads
)
