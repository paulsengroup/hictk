# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT WIN32)
  file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)
endif()

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/15005955/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=4e1ed0871b07c8724d3db635ba45850c1b97a2ee18dcee947d827ca5c5b4c1d2
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# gersemi: on

file(
  ARCHIVE_EXTRACT
  INPUT "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION "${PROJECT_SOURCE_DIR}"
)
