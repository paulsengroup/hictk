# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat REQUIRED)

add_library(hictk_internal_numeric INTERFACE)
add_library(hictk::numeric ALIAS hictk_internal_numeric)
target_sources(
  hictk_internal_numeric
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_numeric
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_system_libraries(hictk_internal_numeric INTERFACE FastFloat::fast_float)
