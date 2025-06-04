# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_tools_definitions INTERFACE)

target_compile_definitions(
  hictk_internal_tools_definitions
  INTERFACE
    "SPDLOG_ACTIVE_LEVEL=$<IF:$<CONFIG:Debug>,SPDLOG_LEVEL_DEBUG,SPDLOG_LEVEL_INFO>"
    SPDLOG_CLOCK_COARSE
    SPDLOG_NO_THREAD_ID
    SPDLOG_DISABLE_DEFAULT_LOGGER
)

add_subdirectory(balance)
add_subdirectory(cli)
add_subdirectory(compressed_io)
add_subdirectory(config)
add_subdirectory(convert)
add_subdirectory(dump)
add_subdirectory(file_attributes_formatting)
add_subdirectory(fix_mcool)
add_subdirectory(load)
add_subdirectory(merge)
add_subdirectory(metadata)
add_subdirectory(rename_chromosomes)
add_subdirectory(telemetry)
add_subdirectory(validate)
add_subdirectory(zoomify)

find_package(CLI11 REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

add_executable(hictk)

target_sources(hictk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

target_include_directories(hictk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_config
    hictk_internal_tools_definitions
    hictk_internal_tools_cli
    hictk_internal_tools_telemetry
    hictk_internal_tools_balance
    hictk_internal_tools_convert
    hictk_internal_tools_dump
    hictk_internal_tools_fix_mcool
    hictk_internal_tools_load
    hictk_internal_tools_merge
    hictk_internal_tools_metadata
    hictk_internal_tools_rename_chromosomes
    hictk_internal_tools_validate
    hictk_internal_tools_zoomify
    hictk::version
)

target_link_system_libraries(
  hictk
  PRIVATE
    CLI11::CLI11
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
  PUBLIC
    Threads::Threads
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    hictk
    PROPERTIES
      LINK_FLAGS_RELEASE
        -s
  )
endif()

include(GNUInstallDirs)
install(
  TARGETS
    hictk
  COMPONENT Runtime
  LIBRARY
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
  PUBLIC_HEADER
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
  PRIVATE_HEADER
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/scratch/"
)

install(FILES "${PROJECT_SOURCE_DIR}/LICENSE" DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/hictk/" COMPONENT Runtime)
