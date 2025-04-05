# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_executable(hictk_pixel_formatting_bench pixel_formatting.cpp)

target_link_libraries(
  hictk_pixel_formatting_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::format
    hictk::pixel
)

target_link_system_libraries(
  hictk_pixel_formatting_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
