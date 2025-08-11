# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(libdeflate REQUIRED QUIET)

add_library(hictk_internal_common_default_delete_libdeflate INTERFACE)
add_library(hictk::default_delete_libdeflate ALIAS hictk_internal_common_default_delete_libdeflate)

target_sources(
  hictk_internal_common_default_delete_libdeflate
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/default_delete_libdeflate.hpp"
)

target_include_directories(
  hictk_internal_common_default_delete_libdeflate
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_system_libraries(
  hictk_internal_common_default_delete_libdeflate
  INTERFACE
    "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)
