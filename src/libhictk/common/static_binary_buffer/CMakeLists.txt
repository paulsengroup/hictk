# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_common_static_binary_buffer INTERFACE)
add_library(hictk::static_binary_buffer ALIAS hictk_internal_common_static_binary_buffer)

target_sources(
  hictk_internal_common_static_binary_buffer
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/static_binary_buffer.hpp"
)

target_include_directories(
  hictk_internal_common_static_binary_buffer
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
