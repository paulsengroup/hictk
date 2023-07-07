# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)

add_executable(hictk_test_compare_coolers)

target_sources(hictk_test_compare_coolers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/hictk_test_compare_coolers.cpp)

target_link_libraries(hictk_test_compare_coolers PRIVATE hictk_project_warnings hictk_project_options hictk::cooler)

target_link_system_libraries(
  hictk_test_compare_coolers
  PRIVATE
  CLI11::CLI11
  std::filesystem)
