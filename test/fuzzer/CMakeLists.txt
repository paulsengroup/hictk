# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

message(STATUS "Building hictk_fuzzer")

add_library(hictk_fuzzer_project_options INTERFACE)
add_library(hictk_fuzzer_project_warnings INTERFACE)
add_library(hictk_fuzzer::project_options ALIAS hictk_fuzzer_project_options)
add_library(hictk_fuzzer::project_warnings ALIAS hictk_fuzzer_project_warnings)
target_link_libraries(hictk_fuzzer_project_options INTERFACE hictk::project_options)
target_link_libraries(hictk_fuzzer_project_warnings INTERFACE hictk::project_warnings)

if(NOT HICTK_WITH_EIGEN OR NOT HICTK_WITH_ARROW)
  message(
    FATAL_ERROR
    "HICTK_ENABLE_FUZZY_TESTING=ON requires HICTK_WITH_EIGEN and HICTK_WITH_ARROW to also be set to ON"
  )
endif()

add_subdirectory(src)
