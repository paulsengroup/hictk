# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15012705/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=d8114081f46eefdf8bea9655e9a6917c28bb46f641244cba2b60d9553d18c3dd
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# gersemi: on

file(
  ARCHIVE_EXTRACT
  INPUT "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION "${PROJECT_SOURCE_DIR}"
)
