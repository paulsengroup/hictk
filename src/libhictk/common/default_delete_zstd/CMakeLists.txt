# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(zstd REQUIRED QUIET)

add_library(hictk_internal_common_default_delete_zstd INTERFACE)
add_library(hictk::common_default_delete_zstd ALIAS hictk_internal_common_default_delete_zstd)

target_sources(
  hictk_internal_common_default_delete_zstd
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_common_default_delete_zstd
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_system_libraries(
  hictk_internal_common_default_delete_zstd
  INTERFACE
    "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)
