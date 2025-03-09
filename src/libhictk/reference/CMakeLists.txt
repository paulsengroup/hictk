# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap REQUIRED)

add_library(hictk_internal_reference INTERFACE)
add_library(hictk::reference ALIAS hictk_internal_reference)

target_sources(
  hictk_internal_reference
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_include_directories(
  hictk_internal_reference
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_reference
  INTERFACE
    hictk::chromosome
    hictk::common
    hictk::format
)

target_link_system_libraries(hictk_internal_reference INTERFACE phmap)
