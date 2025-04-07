# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_executable(hictk_bin_table_coords_to_id_bench bin_table_coords_to_id.cpp)

target_link_libraries(
  hictk_bin_table_coords_to_id_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::bin_table
)

target_link_system_libraries(
  hictk_bin_table_coords_to_id_bench
  PUBLIC
    CLI11::CLI11
    std::filesystem
)
