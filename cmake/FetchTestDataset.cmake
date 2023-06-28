# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/l6rymg9mezixin6/hictk_test_data.tar.xz?dl=1
  EXPECTED_HASH SHA256=13e9a0d293574ad6a2d88f68a1226c5fd69df3dff8abbc5f89c98a1f11421a46
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
