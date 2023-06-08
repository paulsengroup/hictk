# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 CONFIG REQUIRED)
find_package(Filesystem REQUIRED)
find_package(fmt CONFIG REQUIRED)

include(CTest)
include(Catch)

add_executable(coolerpp_test_main)

target_sources(
  coolerpp_test_main
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/attribute_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/balancing_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/bin_table_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/chromosome_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/file_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/dataset_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/index_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/pixel_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/pixel_selector_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/utils_merge_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/utils_equal_test.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/units/variant_buff_test.cpp)

target_include_directories(coolerpp_test_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/units/include/)

target_link_libraries(
  coolerpp_test_main
  PRIVATE coolerpp_project_warnings coolerpp_project_options
  PUBLIC Coolerpp::Coolerpp fmt::fmt)

target_link_system_libraries(
  coolerpp_test_main
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Testing/)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  coolerpp_test_main
  TEST_SPEC
  "[short]"
  TEST_SUFFIX
  " - SHORT"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  OUTPUT_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/Testing/
  EXTRA_ARGS
  --success
  --skip-benchmarks)

catch_discover_tests(
  coolerpp_test_main
  TEST_SPEC
  "[medium]"
  TEST_SUFFIX
  " - MEDIUM"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  EXTRA_ARGS
  --success
  --skip-benchmarks)

catch_discover_tests(
  coolerpp_test_main
  TEST_SPEC
  "[long]"
  TEST_SUFFIX
  " - LONG"
  WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/../
  EXTRA_ARGS
  --success
  --skip-benchmarks)
