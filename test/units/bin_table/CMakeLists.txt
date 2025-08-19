# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)
find_package(Catch2 REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)

include(CTest)
include(Catch)

add_executable(bin_table_tests)

target_sources(
  bin_table_tests
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/bin_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/bin_table_test.cpp"
)

target_link_libraries(
  bin_table_tests
  PRIVATE
    hictk::project_warnings
    hictk::project_options
    hictk::bin_table
    hictk::format
)

target_link_system_libraries(
  bin_table_tests
  PRIVATE
    Catch2::Catch2WithMain
    fmt::fmt-header-only
    std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  bin_table_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
