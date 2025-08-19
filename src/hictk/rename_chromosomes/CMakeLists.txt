# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_rename_chromosomes STATIC)

target_sources(hictk_internal_tools_rename_chromosomes PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/rename_chromosomes.cpp")

target_include_directories(hictk_internal_tools_rename_chromosomes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_rename_chromosomes
  PRIVATE
    hictk::project_options
    hictk::project_warnings
    hictk::cooler
    hictk::hic
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_rename_chromosomes
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
    phmap
)
