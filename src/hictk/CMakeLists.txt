# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(FMT REQUIRED)
find_package(LibArchive REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(readerwriterqueue REQUIRED)
find_package(spdlog REQUIRED)
find_package(tomlplusplus REQUIRED)
find_package(Filesystem REQUIRED)
find_package(Threads REQUIRED)

add_executable(hictk)

target_sources(
  hictk
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_balance.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_convert.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_dump.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_fix_mcool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_load.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_merge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_metadata.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_rename_chromosomes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_validate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_zoomify.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/balance/balance_ice.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/balance/balance_scale.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/balance/balance_vc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/attributes_formatting.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/compressed_reader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/convert/convert.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/convert/cool_to_hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/convert/hic_to_cool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dump/dump.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dump/dump_common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fix_mcool/fix_mcool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load/common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load/hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load/init_bin_table.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load/load.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load/pixel_parser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/merge/cool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/merge/hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/merge/merge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/cool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/mcool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/scool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/metadata/metadata.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rename_chromosomes/rename_chromosomes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validate/cool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validate/hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validate/mcool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validate/scool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/validate/validate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/zoomify/cool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/zoomify/hic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/zoomify/mcool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/zoomify/zoomify.cpp
)

target_compile_definitions(
  hictk
  PRIVATE
    SPDLOG_ACTIVE_LEVEL=$<IF:$<CONFIG:Debug>,SPDLOG_LEVEL_DEBUG,SPDLOG_LEVEL_INFO>
    SPDLOG_CLOCK_COARSE
    SPDLOG_NO_THREAD_ID
    SPDLOG_DISABLE_DEFAULT_LOGGER
)
target_include_directories(hictk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(
  hictk
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk::libhictk
)

target_link_system_libraries(
  hictk
  PRIVATE
  CLI11::CLI11
  fmt::fmt-header-only
  nlohmann_json::nlohmann_json
  readerwriterqueue::readerwriterqueue
  spdlog::spdlog_header_only
  std::filesystem
  tomlplusplus::tomlplusplus
  PUBLIC
  LibArchive::LibArchive
  Threads::Threads
)

include(GNUInstallDirs)
install(
  TARGETS
    hictk
  COMPONENT Runtime
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/
  PRIVATE_HEADER
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/
)

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/hictk/ COMPONENT Runtime)
