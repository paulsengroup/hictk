# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)

add_library(hictk_reference INTERFACE)
add_library(hictk::reference ALIAS hictk_reference)

target_sources(
  hictk_reference
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/reference_impl.hpp)
target_include_directories(hictk_reference INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(
  hictk_reference
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::chromosome
            hictk::common
            hictk::fmt)

target_link_system_libraries(
  hictk_reference
  INTERFACE
  std::filesystem
  _hictk_phmap_tgt)

install(TARGETS hictk_reference FILE_SET HEADERS EXCLUDE_FROM_ALL)
