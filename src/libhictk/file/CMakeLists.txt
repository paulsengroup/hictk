# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT CONFIG REQUIRED)

add_library(file INTERFACE)
add_library(hictk::file ALIAS file)

target_sources(
  file
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(file INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                          "$<INSTALL_INTERFACE:include>")
target_link_libraries(
  file
  INTERFACE hictk::balancing
            hictk::cooler
            hictk::hic
            hictk::reference)

target_link_system_libraries(file INTERFACE fmt::fmt)
