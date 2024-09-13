# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/13760444/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=621923992c601e2d9c72975b10c56d2ea3fb507e9c0c69150165c7ea16883051
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
