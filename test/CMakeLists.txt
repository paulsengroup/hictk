# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(units)

if(HICTK_ENABLE_FUZZY_TESTING)
  add_subdirectory(fuzzer)
endif()

if(HICTK_DOWNLOAD_TEST_DATASET)
  include("${PROJECT_SOURCE_DIR}/cmake/FetchTestDataset.cmake")
endif()
