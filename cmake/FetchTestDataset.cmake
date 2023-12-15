# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/10289491/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=5e69dceb8789d923a38aed7add8fc18abfdfe531aea6effcdb7efe3c9bcf5246
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz")
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
