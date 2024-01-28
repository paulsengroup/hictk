# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(units)
add_subdirectory(utils)

option(HICTK_DOWNLOAD_TEST_DATASET "Download datasets required by unit and integration tests" ON)

if(HICTK_DOWNLOAD_TEST_DATASET)
  message(STATUS "Downloading test dataset...")
  include(${PROJECT_SOURCE_DIR}/cmake/FetchTestDataset.cmake)
endif()
