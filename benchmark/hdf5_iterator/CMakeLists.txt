# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(hictk_hdf5_iterator_bench hdf5_iterator.cpp)

target_link_libraries(
  hictk_hdf5_iterator_bench
  PRIVATE
    hictk_project_warnings
    hictk_project_options
  PUBLIC
    hictk::cooler
)

target_link_system_libraries(
  hictk_hdf5_iterator_bench
  PUBLIC
  CLI11::CLI11
  std::filesystem
)
