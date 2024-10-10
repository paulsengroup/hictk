# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(hictk_zoomify_bench zoomify.cpp)

target_link_libraries(
  hictk_zoomify_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::transformers
)

target_link_system_libraries(
  hictk_zoomify_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
