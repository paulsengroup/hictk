# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT CONFIG REQUIRED)

add_library(formatting INTERFACE)
add_library(hictk::format ALIAS formatting)

target_sources(
  formatting
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(formatting INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                                $<INSTALL_INTERFACE:include>)

target_link_libraries(
  formatting
  INTERFACE hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::genomic_interval
            hictk::pixel)

target_link_system_libraries(formatting INTERFACE fmt::fmt)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
