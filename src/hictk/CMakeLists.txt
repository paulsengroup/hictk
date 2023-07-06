# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Boost COMPONENTS headers filesystem)
if(NOT boost_FOUND)
  find_package(Boost CONFIG REQUIRED COMPONENTS headers filesystem)
endif()

find_package(Filesystem REQUIRED)
find_package(Threads REQUIRED)

add_executable(hictk_tools)
set_target_properties(hictk_tools PROPERTIES OUTPUT_NAME hictk)

target_sources(
  hictk_tools
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/cli/cli.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/convert.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/cool_to_hic.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/convert/hic_to_cool.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dump/dump.cpp)

target_include_directories(hictk_tools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
  hictk_tools
  PRIVATE hictk_project_options
          hictk_project_warnings
          hictk::cooler
          hictk::hic)

target_link_system_libraries(
  hictk_tools
  PRIVATE
  Boost::headers
  CLI11::CLI11
  readerwriterqueue
  std::filesystem
  PUBLIC
  Boost::filesystem
  fmt::fmt
  spdlog::spdlog
  Threads::Threads)

include(GNUInstallDirs)
install(
  TARGETS hictk_tools
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/hictk/)
