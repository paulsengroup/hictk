# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(nlohmann_json REQUIRED QUIET)
find_package(tomlplusplus REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_file_attributes_formatting STATIC)

target_sources(
  hictk_internal_tools_file_attributes_formatting
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/file_attributes_formatting.cpp"
)

target_include_directories(hictk_internal_tools_file_attributes_formatting PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_file_attributes_formatting
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk::common
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_file_attributes_formatting
  PRIVATE
  fmt::fmt-header-only
    std::filesystem
  PUBLIC
  nlohmann_json::nlohmann_json
  tomlplusplus::tomlplusplus
)
