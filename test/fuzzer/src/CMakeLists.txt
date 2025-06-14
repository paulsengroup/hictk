# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(PYBIND11_FINDPYTHON ON)
find_package(Boost REQUIRED QUIET)
find_package(CLI11 REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(pybind11 REQUIRED)
find_package(spdlog REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

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
target_compile_definitions(
  hictk_fuzzer
  PRIVATE
    SPDLOG_ACTIVE_LEVEL=$<IF:$<CONFIG:Debug>,SPDLOG_LEVEL_DEBUG,SPDLOG_LEVEL_INFO>
    SPDLOG_CLOCK_COARSE
    SPDLOG_NO_THREAD_ID
    SPDLOG_DISABLE_DEFAULT_LOGGER
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
    "$<TARGET_NAME_IF_EXISTS:Boost::process>"
    CLI11::CLI11
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
  PRIVATE
    pybind11::embed
    Threads::Threads
)
