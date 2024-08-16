# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(hictk_cooler_creation_bench cooler_creation.cpp)

target_link_libraries(
  hictk_cooler_creation_bench
  PRIVATE hictk_project_warnings hictk_project_options
  PUBLIC hictk::cooler)

target_link_system_libraries(
  hictk_cooler_creation_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem)
