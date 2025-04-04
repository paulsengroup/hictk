# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

add_library(hictk_internal_tools_cli STATIC)

target_sources(
  hictk_internal_tools_cli
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_balance.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_convert.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_dump.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_fix_mcool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_load.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_merge.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_metadata.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_rename_chromosomes.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_validate.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cli_zoomify.cpp"
)

target_include_directories(hictk_internal_tools_cli PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_cli
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk_internal_tools_telemetry
    hictk::common
    hictk::cooler
    hictk::file
    hictk::hic
    hictk::tmpdir
  PUBLIC
    hictk_internal_tools_config
    hictk::common_genomic_units # TODO fixme: this should be private
)

target_link_system_libraries(
  hictk_internal_tools_cli
  PRIVATE
  CLI11::CLI11
  fmt::fmt-header-only
  spdlog::spdlog_header_only
  std::filesystem
  phmap
  Threads::Threads
)
