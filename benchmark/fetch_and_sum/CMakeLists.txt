# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_executable(hictk_fetch_and_sum_bench fetch_and_sum.cpp)

target_link_libraries(
  hictk_fetch_and_sum_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
    hictk::hic
)

target_link_system_libraries(
  hictk_fetch_and_sum_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
