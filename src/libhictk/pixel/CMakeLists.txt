# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_pixel INTERFACE)
add_library(hictk::pixel ALIAS hictk_pixel)

target_sources(
  hictk_pixel
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/pixel_impl.hpp)
target_include_directories(hictk_pixel INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(
  hictk_pixel
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::fmt)

install(TARGETS hictk_pixel FILE_SET HEADERS EXCLUDE_FROM_ALL)
