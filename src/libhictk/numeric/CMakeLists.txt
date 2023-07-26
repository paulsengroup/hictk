# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FastFloat CONFIG REQUIRED)

add_library(hictk_numeric INTERFACE)
add_library(hictk::numeric ALIAS hictk_numeric)
target_sources(
  hictk_numeric
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/numeric_utils_impl.hpp)

target_include_directories(hictk_numeric INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(hictk_numeric INTERFACE hictk_project_warnings hictk_project_options)

target_link_system_libraries(hictk_numeric INTERFACE FastFloat::fast_float)

install(TARGETS hictk_numeric FILE_SET HEADERS EXCLUDE_FROM_ALL)
