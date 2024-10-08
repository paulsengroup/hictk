# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(Catch2 REQUIRED)
find_package(fmt REQUIRED)

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
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::bin_table
    hictk::format
)

target_link_system_libraries(
  bin_table_tests
  PUBLIC
  Catch2::Catch2WithMain
  fmt::fmt-header-only
  std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  bin_table_tests
  TEST_SPEC
  "[short]"
  TEST_SUFFIX
  " - SHORT"
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --success
  --skip-benchmarks
)

catch_discover_tests(
  bin_table_tests
  TEST_SPEC
  "[medium]"
  TEST_SUFFIX
  " - MEDIUM"
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  EXTRA_ARGS
  --success
  --skip-benchmarks
)

catch_discover_tests(
  bin_table_tests
  TEST_SPEC
  "[long]"
  TEST_SUFFIX
  " - LONG"
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  EXTRA_ARGS
  --success
  --skip-benchmarks
)
