# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/10790574/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=9f6f1ee8c8bfe3c904f55f58cbdf51730d6f3fa57be31ddb8c51e81ee0872fb5
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
