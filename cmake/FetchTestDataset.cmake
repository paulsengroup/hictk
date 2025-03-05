# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/test/data/")
file(LOCK "${PROJECT_SOURCE_DIR}/test/data/" DIRECTORY GUARD FILE)

# gersemi: off
file(
  DOWNLOAD https://zenodo.org/records/13851354/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=55071172638948112a69a43ebfd07b0b9b830cc5e1bfef87b05b586d228ab1bd
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# gersemi: on

file(
  ARCHIVE_EXTRACT
  INPUT "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION "${PROJECT_SOURCE_DIR}"
)
