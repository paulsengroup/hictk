# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 REQUIRED)
find_package(Filesystem REQUIRED)
find_package(FMT REQUIRED)
include(CTest)
include(Catch)

add_executable(hictk_formatting_bench)

target_sources(hictk_formatting_bench PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/formatting.cpp")

target_link_libraries(
  hictk_formatting_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::bin_table
    hictk::chromosome
    hictk::format
    hictk::pixel
    hictk::reference
)

target_link_system_libraries(
  hictk_formatting_bench
  PRIVATE fmt::fmt-header-only
  PUBLIC Catch2::Catch2WithMain std::filesystem
)

catch_discover_tests(
  hictk_formatting_bench
  TEST_SPEC
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Benchmarking/"
  EXTRA_ARGS
  --success
)
