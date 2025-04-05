# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_executable(hictk_hic_traversal_bench hic_traversal.cpp)

target_link_libraries(
  hictk_hic_traversal_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::hic
)

target_link_system_libraries(
  hictk_hic_traversal_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
