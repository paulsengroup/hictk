# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(hictk_bin_table_id_to_coords_bench bin_table_id_to_coords.cpp)

target_link_libraries(
  hictk_bin_table_id_to_coords_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::bin_table
)

target_link_system_libraries(
  hictk_bin_table_id_to_coords_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
