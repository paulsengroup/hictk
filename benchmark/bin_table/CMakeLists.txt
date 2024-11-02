# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 REQUIRED)
find_package(Filesystem REQUIRED)
find_package(FMT REQUIRED)
include(CTest)
include(Catch)

add_executable(hictk_bin_table_bench)
add_dependencies(hictk_benchmarks hictk_bin_table_bench)

target_sources(
  hictk_bin_table_bench
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/bin_ids_to_genomic_coords.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/genomic_coords_to_bin_ids.cpp"
)

target_link_libraries(
  hictk_bin_table_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::bin_table
)

target_link_system_libraries(
  hictk_bin_table_bench
  PRIVATE fmt::fmt-header-only
  PUBLIC Catch2::Catch2WithMain std::filesystem
)

catch_discover_tests(
  hictk_bin_table_bench
  TEST_SPEC
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Benchmarking/"
  EXTRA_ARGS
  --success
)
