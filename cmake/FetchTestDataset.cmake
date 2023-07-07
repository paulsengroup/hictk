# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/8124768/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=e9988753271d966731d0af1c95f5a043227aea7b36643bbe91730483b7056cab
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
