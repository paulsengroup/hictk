# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT CONFIG REQUIRED)

add_library(hictk_fmt INTERFACE)
add_library(hictk::fmt ALIAS hictk_fmt)

target_sources(
  hictk_fmt
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(hictk_fmt INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                               $<INSTALL_INTERFACE:include>)

target_link_libraries(
  hictk_fmt
  INTERFACE hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::genomic_interval
            hictk::pixel)

target_link_system_libraries(hictk_fmt INTERFACE fmt::fmt)

include(GNUInstallDirs)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
