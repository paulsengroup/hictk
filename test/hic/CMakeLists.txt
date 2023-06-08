# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

find_package(
    Catch2
    3
    CONFIG
    QUIET)

include(CTest)
include(Catch)

add_executable(hicxx_test_main)

target_sources(
    hicxx_test_main
    PRIVATE filestream_local_test.cpp
            filestream_remote_test.cpp
            hic_file_test.cpp
            hic_file_stream_test.cpp
            matrix_zoom_data_test.cpp)

target_link_libraries(
    hicxx_test_main
    PRIVATE hicxx_project_options
            hicxx_project_warnings
            hicxx
            Catch2::Catch2WithMain)

catch_discover_tests(hicxx_test_main WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../)
