# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FetchContent_Declare(
  _hictk_catch2
  URL ${PROJECT_SOURCE_DIR}/external/catch2-v3.3.2.tar.xz
  URL_HASH SHA256=05ac83d7c65d6ee4a499b8c98b32ef2cb30280d7691037981feffe498992804f
  SYSTEM)

FetchContent_MakeAvailable(_hictk_catch2)

list(APPEND CMAKE_MODULE_PATH ${_hictk_catch2_SOURCE_DIR}/extras)

add_subdirectory(units)

option(HICTK_DOWNLOAD_TEST_DATASET "Download datasets required by unit and integration tests" ON)

if(HICTK_DOWNLOAD_TEST_DATASET)
  message(STATUS "Downloading test dataset...")
  include(${PROJECT_SOURCE_DIR}/cmake/FetchTestDataset.cmake)
endif()
