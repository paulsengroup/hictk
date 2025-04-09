# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

set(TEST_DATASET_TAR "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst")

message(STATUS "Fetching the test dataset")

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15182021/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=7734a78307885e549ca64828a497a2c524e47acb92453e42532409ec7fb4e98a
  "${TEST_DATASET_TAR}"
)
# gersemi: on

message(STATUS "Fetching the test dataset - done")

message(STATUS "Extracting the test dataset")

file(ARCHIVE_EXTRACT INPUT "${TEST_DATASET_TAR}" DESTINATION "${PROJECT_SOURCE_DIR}")

message(STATUS "Extracting the test dataset - done")
message(STATUS "Test datasets can be found under \"${PROJECT_SOURCE_DIR}/test/data/\"")
