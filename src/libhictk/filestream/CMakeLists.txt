# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_filestream STATIC)
add_library(hictk::filestream ALIAS hictk_internal_filestream)

target_sources(
  hictk_internal_filestream
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/filestream.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/filestream_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/filestream.cpp"
)

target_include_directories(
  hictk_internal_filestream
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_filestream
  PUBLIC
    hictk::common
    hictk::type_traits
)

target_link_system_libraries(hictk_internal_filestream PRIVATE fmt::fmt-header-only)
