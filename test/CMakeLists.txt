# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FetchContent_Declare(
  _hictk_catch2
  URL ${PROJECT_SOURCE_DIR}/external/catch2-v3.5.2.tar.xz
  URL_HASH SHA256=b680c9665cb6527e012248f1bfc4b50c90e4a03f88874927bc9e06d3cb14bf7b
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_MakeAvailable(_hictk_catch2)
target_disable_clang_tidy(Catch2)
target_disable_clang_tidy(Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${_hictk_catch2_SOURCE_DIR}/extras)

add_subdirectory(units)
add_subdirectory(utils)

option(HICTK_DOWNLOAD_TEST_DATASET "Download datasets required by unit and integration tests" ON)

if(HICTK_DOWNLOAD_TEST_DATASET)
  message(STATUS "Downloading test dataset...")
  include(${PROJECT_SOURCE_DIR}/cmake/FetchTestDataset.cmake)
endif()
