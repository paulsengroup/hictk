# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictk_hic_tests)

target_sources(
  hictk_hic_tests
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/file_reader_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_writer_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_zoomify_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_file_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_merge_test.cpp"
)

target_link_libraries(
  hictk_hic_tests
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::hic
    hictk::tmpdir
)

target_link_system_libraries(
  hictk_hic_tests
  PUBLIC
    Catch2::Catch2WithMain
    std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictk_hic_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
