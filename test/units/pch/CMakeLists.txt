# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED)

if(HICTK_USE_PCH)
  add_library(hictk_test_pch)
  target_sources(hictk_test_pch PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/lib.cpp")
  target_precompile_headers(hictk_test_pch PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/pch.hpp")
  target_link_libraries(
    hictk_test_pch
    PRIVATE hictk_project_options hictk_project_warnings
    PUBLIC hictk::libhictk
           Catch2::Catch2WithMain
           fmt::fmt-header-only)
endif()
