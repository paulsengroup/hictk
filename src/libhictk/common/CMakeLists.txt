# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)
find_package(libdeflate REQUIRED)
find_package(spdlog REQUIRED)
find_package(zstd REQUIRED)

add_library(hictk_internal_common INTERFACE)
add_library(hictk::common ALIAS hictk_internal_common)
add_dependencies(hictk_internal_common _hictk_check_git)

target_sources(
  hictk_internal_common
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_common
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_common
  INTERFACE
    fmt::fmt-header-only
    "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
    spdlog::spdlog_header_only
    "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)
