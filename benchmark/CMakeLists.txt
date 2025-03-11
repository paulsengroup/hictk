# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Catch2 REQUIRED)
find_package(Filesystem REQUIRED)
find_package(FMT REQUIRED)
include(CTest)
include(Catch)

add_executable(hictk_benchmark)

target_sources(hictk_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

add_subdirectory(bin_table)
add_subdirectory(bins)
add_subdirectory(cooler_creation)
add_subdirectory(cooler_traversal)
add_subdirectory(formatting)
add_subdirectory(hdf5_iterator)
add_subdirectory(hic_creation)
add_subdirectory(hic_traversal)
add_subdirectory(interaction_fetching)
add_subdirectory(pixel_merger)
add_subdirectory(pixels)
add_subdirectory(reference)
add_subdirectory(zoomify)

target_include_directories(hictk_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/")

target_link_libraries(
  hictk_benchmark
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::file
    hictk::hic
)

target_link_system_libraries(
        hictk_benchmark
        PRIVATE fmt::fmt-header-only std::filesystem
        PUBLIC Catch2::Catch2
)

catch_discover_tests(
        hictk_benchmark
        TEST_SPEC
        WORKING_DIRECTORY
        "${PROJECT_SOURCE_DIR}"
        OUTPUT_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/Benchmarking/"
        EXTRA_ARGS
        --success
)
