# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/records/13825195/files/hictk_test_data.tar.zst?download=1
  EXPECTED_HASH SHA256=e5a93842a3d068d752f28da8cc0b9ba90064ad81ccfaa811fd68e96c5af221dc
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  "${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.zst"
  DESTINATION
  "${PROJECT_SOURCE_DIR}")
