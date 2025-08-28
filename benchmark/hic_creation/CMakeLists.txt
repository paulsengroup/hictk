# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_executable(hictk_hic_creation_bench hic_creation.cpp)

target_link_libraries(
  hictk_hic_creation_bench
  PRIVATE
    hictk::project_warnings
    hictk::project_options
  PUBLIC
    hictk::hic
    hictk::hic
)

target_link_system_libraries(
  hictk_hic_creation_bench
  PUBLIC
    CLI11::CLI11
    std::filesystem
)
