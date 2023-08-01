# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)

include(CTest)
include(Catch)

add_executable(hictk_variant_tests)

target_sources(hictk_variant_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/variant_test.cpp")

target_link_libraries(
  hictk_variant_tests
  PRIVATE hictk_project_warnings hictk_project_options
  PUBLIC hictk::variant)

target_link_system_libraries(
  hictk_variant_tests
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  hictk_variant_tests
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
  --skip-benchmarks)

catch_discover_tests(
  hictk_variant_tests
  TEST_SPEC
  "[medium]"
  TEST_SUFFIX
  " - MEDIUM"
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  EXTRA_ARGS
  --success
  --skip-benchmarks)

catch_discover_tests(
  hictk_variant_tests
  TEST_SPEC
  "[long]"
  TEST_SUFFIX
  " - LONG"
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  EXTRA_ARGS
  --success
  --skip-benchmarks)
