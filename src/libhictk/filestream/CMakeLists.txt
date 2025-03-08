# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap REQUIRED)
find_package(spdlog REQUIRED)

add_library(hictk_filestream INTERFACE)
add_library(hictk::filestream ALIAS hictk_filestream)

target_sources(
  hictk_filestream
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_filestream
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(hictk_filestream INTERFACE hictk::common)
