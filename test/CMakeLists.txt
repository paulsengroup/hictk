# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FetchContent_Declare(
  _hictk_catch2
  URL ${PROJECT_SOURCE_DIR}/external/catch2-v3.4.0.tar.xz
  URL_HASH SHA256=757d1b29b836e2c816784e062785e57ad852317de6242465129d4acb65b58b4d
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
