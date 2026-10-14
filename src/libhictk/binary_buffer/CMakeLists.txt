# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_binary_buffer STATIC)
add_library(hictk::binary_buffer ALIAS hictk_internal_binary_buffer)

target_sources(
  hictk_internal_binary_buffer
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/binary_buffer.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/binary_buffer_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/binary_buffer.cpp"
)

target_include_directories(
  hictk_internal_binary_buffer
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
