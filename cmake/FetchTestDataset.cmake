# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/8121687/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=4ac996d0169410036aa1c8a0aeb49da9c44a44cf47b255eb589650adc5ea3468
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
