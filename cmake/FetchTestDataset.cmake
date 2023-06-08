# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/6svpdoeh29qyebt/seqtk_test_data.tar.xz?dl=1
  EXPECTED_HASH SHA512=ac65bd5afc6ef76ca8ef1ce8c4cb1028c769cf47e6f78d715c51bd5bc015520f3cb7a0c1a37f6b3ddc85cc4e921ff3d64c03de68130a79e85ff742cc33d8787f
  ${PROJECT_SOURCE_DIR}/test/data/seqtk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/seqtk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
