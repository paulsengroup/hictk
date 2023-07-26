# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_balancing INTERFACE)
add_library(hictk::balancing ALIAS hictk_balancing)

target_sources(
  hictk_balancing
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/weights_impl.hpp)

target_include_directories(hictk_balancing INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(
  hictk_balancing
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::common
            hictk::pixel)

target_link_system_libraries(hictk_balancing INTERFACE _hictk_phmap_tgt)

install(TARGETS hictk_balancing FILE_SET HEADERS EXCLUDE_FROM_ALL)
