# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15068509/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=85960cc4088c8bdf6122f9c6660cefcf189e0fb050e30090988ae719c9ddb19f
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# gersemi: on

file(
  ARCHIVE_EXTRACT
  INPUT "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION "${PROJECT_SOURCE_DIR}"
)
