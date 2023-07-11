# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://zenodo.org/record/8136686/files/hictk_test_data.tar.xz?download=1
  EXPECTED_HASH SHA256=ce4b9c08a2715f1a43be7cc5f3d44e066508088fd64363d040222e316f030f1d
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/hictk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
