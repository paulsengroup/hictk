# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(PYBIND11_FINDPYTHON ON)
find_package(Boost REQUIRED)
find_package(CLI11 REQUIRED)
find_package(Filesystem REQUIRED)
find_package(FMT REQUIRED)
find_package(pybind11 REQUIRED)
find_package(spdlog REQUIRED)
find_package(Threads REQUIRED)

add_executable(hictk_fuzzer)

target_sources(
  hictk_fuzzer
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cli.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/fuzzer_scheduler.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/fuzzer_worker.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)

target_include_directories(hictk_fuzzer PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(
  hictk_fuzzer
  PROPERTIES
    CXX_VISIBILITY_PRESET
      hidden
)
target_link_libraries(
  hictk_fuzzer
  PRIVATE
    hictk_fuzzer_project_options
    hictk_fuzzer_project_warnings
    hictk::libhictk
)

target_link_system_libraries(
  hictk_fuzzer
  PRIVATE
  Boost::headers
  CLI11::CLI11
  fmt::fmt-header-only
  spdlog::spdlog_header_only
  std::filesystem
  PRIVATE
  pybind11::embed
  Threads::Threads
)
