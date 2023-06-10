# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        _hictk_catch2
        URL ${PROJECT_SOURCE_DIR}/external/catch2-v3.3.2.tar.xz
        URL_HASH SHA512=a85f6dc72da438034ab46a2a026b9cf703ece5bd4ea62749aa1d9dce71f01cafc04aafbec7aa2e56baf76f9fdf344ad5dda2e8fc42d4ae87e2cffe8fda6bb6ae
)
FetchContent_MakeAvailable(_hictk_catch2)

add_subdirectory(units)

option(HICTK_DOWNLOAD_TEST_DATASET "Download datasets required by unit and integration tests" ON)

if(HICTK_DOWNLOAD_TEST_DATASET)
  message("-- Downloading test dataset...")
  include(${PROJECT_SOURCE_DIR}/cmake/FetchTestDataset.cmake)
endif()
