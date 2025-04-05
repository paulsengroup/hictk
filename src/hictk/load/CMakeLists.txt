# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(bshoshany-thread-pool REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(readerwriterqueue REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

add_library(hictk_internal_tools_load STATIC)

target_sources(
  hictk_internal_tools_load
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/common.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/init_bin_table.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/load.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_parser.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/common.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cool.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/init_bin_table.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/load.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pairs_aggregator.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_parser.hpp"
)

target_include_directories(hictk_internal_tools_load PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_load
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_compressed_io
    hictk_internal_tools_definitions
    hictk::bin_table
    hictk::common
    hictk::cooler
    hictk::genomic_interval
    hictk::hic
    hictk::pixel
    hictk::reference
    hictk::string
    hictk::tmpdir
    hictk::type_traits
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_load
  PRIVATE
  bshoshany-thread-pool::bshoshany-thread-pool
  fmt::fmt-header-only
  readerwriterqueue::readerwriterqueue
  spdlog::spdlog_header_only
  std::filesystem
  phmap
  PUBLIC
  Threads::Threads
)
