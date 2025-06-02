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
  DOWNLOAD https://zenodo.org/records/15579549/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=af98f7997852c4a896fb19f1084dd482b0136b7e3ced15b02bcf71ed96fb5573
  "${TEST_DATASET_TAR}"
)
# gersemi: on

message(STATUS "Fetching the test dataset - done")

message(STATUS "Extracting the test dataset")

file(ARCHIVE_EXTRACT INPUT "${TEST_DATASET_TAR}" DESTINATION "${PROJECT_SOURCE_DIR}")

message(STATUS "Extracting the test dataset - done")
message(STATUS "Test datasets can be found under \"${PROJECT_SOURCE_DIR}/test/data/\"")
