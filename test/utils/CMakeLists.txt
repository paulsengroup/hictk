# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(FMT REQUIRED)
find_package(Filesystem REQUIRED)

add_executable(hictk_test_dump_weights)

target_sources(hictk_test_dump_weights PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/hictk_test_dump_weights.cpp)

target_link_libraries(hictk_test_dump_weights PRIVATE hictk_project_warnings hictk_project_options hictk::libhictk)

target_link_system_libraries(
  hictk_test_dump_weights
  PRIVATE
  CLI11::CLI11
  fmt::fmt-header-only
  std::filesystem)
