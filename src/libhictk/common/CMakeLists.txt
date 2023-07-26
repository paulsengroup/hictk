# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_common INTERFACE)
add_library(hictk::common ALIAS hictk_common)

target_sources(
  hictk_common
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(hictk_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(hictk_common INTERFACE hictk_project_warnings hictk_project_options)

install(TARGETS hictk_common FILE_SET HEADERS EXCLUDE_FROM_ALL)
