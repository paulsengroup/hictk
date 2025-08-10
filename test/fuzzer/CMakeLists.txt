# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

message(STATUS "Building hictk_fuzzer")

add_library(hictk_fuzzer_project_options INTERFACE)
add_library(hictk_fuzzer_project_warnings INTERFACE)
target_link_libraries(hictk_fuzzer_project_options INTERFACE hictk_project_options)
target_link_libraries(hictk_fuzzer_project_warnings INTERFACE hictk_project_warnings)

target_compile_features(hictk_fuzzer_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

if(NOT HICTK_WITH_EIGEN OR NOT HICTK_WITH_ARROW)
  message(
    FATAL_ERROR
    "HICTK_ENABLE_FUZZY_TESTING=ON requires HICTK_WITH_EIGEN and HICTK_WITH_ARROW to also be set to ON"
  )
endif()

add_subdirectory(src)
