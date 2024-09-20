# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/13820402/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=c9f5b2ca7bd635c8e761f21c0215704a4ac76e869307750faaa01a5eec128ef0
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
