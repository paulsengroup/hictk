# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictk_cooler_tests)

target_sources(
  hictk_cooler_tests
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/attribute_numeric_conversions_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/attribute_read_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/attribute_write_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_accessors_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_cow_chunk_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_attributes_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_iteration_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_large_rw_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_read_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dataset_write_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_accessors_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_attributes_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_bins_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_chromosomes_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_ctors_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_pixels_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_sentinel_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_weights_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/index_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/multires_cooler_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/singlecell_cooler_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_1d_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_2d_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_balancing_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_copy_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_equal_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_merge_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_rename_chromosomes_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils_validation_test.cpp"
)

target_link_libraries(
  hictk_cooler_tests
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::tmpdir
)

target_link_system_libraries(
  hictk_cooler_tests
  PUBLIC
    Catch2::Catch2WithMain
    std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictk_cooler_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
