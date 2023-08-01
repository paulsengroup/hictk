# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/8204309/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=f9d0ccc8f8c7bade7f9049ea0abb6b26c061e6fed9e1e7bb8f2da0e1a402cfa4
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
