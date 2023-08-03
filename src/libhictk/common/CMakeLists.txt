# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)

add_library(common INTERFACE)
add_library(hictk::common ALIAS common)

target_sources(
  common
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(common INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                            "$<INSTALL_INTERFACE:include>")

target_link_libraries(common INTERFACE fmt::fmt-header-only)
