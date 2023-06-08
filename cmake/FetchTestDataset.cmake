# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/6svpdoeh29qyebt/seqtk_test_data.tar.xz?dl=1
  EXPECTED_HASH SHA512=c55a36a1d1bdc319dba69839c9872f7745d5a71214973c00732dbbe01ea38b0e1083d61d2b6a76b4d2992ec3cd151bf190a75f81b1a543a1a750f9d0fb84a35a
  ${PROJECT_SOURCE_DIR}/test/data/seqtk_test_data.tar.xz)
# cmake-format: on

file(
  ARCHIVE_EXTRACT
  INPUT
  ${PROJECT_SOURCE_DIR}/test/data/seqtk_test_data.tar.xz
  DESTINATION
  ${PROJECT_SOURCE_DIR})
