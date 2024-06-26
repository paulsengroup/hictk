# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(CLI11 REQUIRED)
find_package(FMT REQUIRED)
find_package(readerwriterqueue REQUIRED)
find_package(spdlog REQUIRED)
find_package(Filesystem REQUIRED)
find_package(Threads REQUIRED)

add_executable(hictk)

target_sources(
  hictk
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_balance.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_convert.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_dump.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_fix_mcool.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_load.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_merge.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_rename_chromosomes.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_validate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli_zoomify.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/balance/balance.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/convert.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/cool_to_hic.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/hic_to_cool.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dump/dump.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dump/dump_common.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/fix_mcool/fix_mcool.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/load/load.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/merge/merge.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/rename_chromosomes/rename_chromosomes.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/validate/validate.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/zoomify/zoomify.cpp)

target_include_directories(hictk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(hictk PRIVATE hictk_project_options hictk_project_warnings hictk::libhictk)

target_link_system_libraries(
  hictk
  PRIVATE
  CLI11::CLI11
  readerwriterqueue::readerwriterqueue
  std::filesystem
  PUBLIC
  fmt::fmt-header-only
  spdlog::spdlog_header_only
  Threads::Threads)

include(GNUInstallDirs)
install(
  TARGETS hictk
  COMPONENT Runtime
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/
  PRIVATE_HEADER DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scratch/)

install(
  FILES ${PROJECT_SOURCE_DIR}/LICENSE
  DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/hictk/
  COMPONENT Runtime)
