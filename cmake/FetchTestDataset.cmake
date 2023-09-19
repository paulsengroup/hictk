# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/l6rymg9mezixin6/hictk_test_data.tar.xz?dl=1
  EXPECTED_HASH SHA256=a97c3a66d25c7441154ef15c9b747e69ac1b6a5810a478a67139565ee3ea999c
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
