# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 REQUIRED)
find_package(Filesystem REQUIRED)
find_package(FMT REQUIRED)
include(CTest)
include(Catch)

add_executable(hictk_interaction_fetching_bench)
add_dependencies(hictk_benchmarks hictk_interaction_fetching_bench)

target_sources(
  hictk_interaction_fetching_bench
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler_trans_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_trans_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_trans_queries.cpp"
)

target_link_libraries(
  hictk_interaction_fetching_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::file
    hictk::hic
)

target_link_system_libraries(
  hictk_interaction_fetching_bench
  PRIVATE fmt::fmt-header-only
  PUBLIC Catch2::Catch2WithMain std::filesystem
)

catch_discover_tests(
  hictk_interaction_fetching_bench
  TEST_SPEC
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Benchmarking/"
  EXTRA_ARGS
  --success
)
