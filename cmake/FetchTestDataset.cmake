# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/13798772/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=0a92335d6adf2991609dea005c77991f4a4c8df0405a6c425eb7176945b7443e
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
