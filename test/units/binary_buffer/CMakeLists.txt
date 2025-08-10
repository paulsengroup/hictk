# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictk_binary_buffer_tests)

target_sources(hictk_binary_buffer_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/binary_buffer_test.cpp")

target_link_libraries(
  hictk_binary_buffer_tests
  PRIVATE
    hictk::project_warnings
    hictk::project_options
    hictk::binary_buffer
)

target_link_system_libraries(
  hictk_binary_buffer_tests
  PRIVATE
    Catch2::Catch2WithMain
    std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictk_binary_buffer_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
