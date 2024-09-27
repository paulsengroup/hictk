# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/13849053/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=4a5a871421a981bc11fb4e4c6b9877e281dd721ef56fc35ff8cb940e81e301a0
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
