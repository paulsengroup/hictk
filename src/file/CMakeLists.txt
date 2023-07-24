# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)

add_library(file INTERFACE)
add_library(hictk::file ALIAS file)

target_sources(
  file
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/file_impl.hpp)

target_include_directories(file INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
  file
  INTERFACE hictk_project_options
            hictk_project_warnings
            hictk::balancing
            hictk::cooler
            hictk::hic
            hictk::reference)

target_link_system_libraries(
  file
  INTERFACE
  _hictk_fmt_tgt
  std::filesystem)
