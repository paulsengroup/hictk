# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

find_package(Catch2 REQUIRED QUIET)
include(CTest)
include(Catch)

add_executable(hictk_transformers_tests)

target_sources(
  hictk_transformers_tests
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_coarsen_pixels_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_diagonal_band_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_join_genomic_coords_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_pixel_merger_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_stats_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_to_dataframe_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_to_dense_matrix_test.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/transformers_to_sparse_matrix_test.cpp"
)

target_link_libraries(
  hictk_transformers_tests
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::tmpdir
    hictk::transformers
)

target_link_system_libraries(
  hictk_transformers_tests
  PUBLIC
  Catch2::Catch2WithMain
  std::filesystem
)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Testing/")

catch_discover_tests(
  hictk_transformers_tests
  WORKING_DIRECTORY
  "${PROJECT_SOURCE_DIR}"
  OUTPUT_DIR
  "${CMAKE_CURRENT_BINARY_DIR}/Testing/"
  EXTRA_ARGS
  --skip-benchmarks
)
