# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictk_balancing_tests)

target_sources(
  hictk_balancing_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/balancing_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/matrix_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/weights_test.cpp
)

target_link_libraries(
  hictk_balancing_tests
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::balancing
    hictk::file
    hictk::tmpdir
)

target_link_system_libraries(
  hictk_balancing_tests
  PUBLIC
    Catch2::Catch2WithMain
    std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictk_balancing_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
