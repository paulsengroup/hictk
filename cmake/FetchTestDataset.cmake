# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/l6rymg9mezixin6/hictk_test_data.tar.xz?dl=1
  EXPECTED_HASH SHA256=ff761ffb1de9cf8dd1437573bac7c0242ef3a75c5621646354e203a9cf9c8397
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
