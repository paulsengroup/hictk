# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(pixel INTERFACE)
add_library(hictk::pixel ALIAS pixel)

target_sources(
  pixel
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(pixel INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                           "$<INSTALL_INTERFACE:include>")
target_link_libraries(
  pixel
  INTERFACE hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::format)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
