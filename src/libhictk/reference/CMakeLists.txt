# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(phmap CONFIG REQUIRED)

add_library(reference INTERFACE)
add_library(hictk::reference ALIAS reference)

target_sources(
  reference
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(reference INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                               $<INSTALL_INTERFACE:include>)
target_link_libraries(reference INTERFACE hictk::chromosome hictk::common hictk::format)

target_link_system_libraries(
  reference
  INTERFACE
  std::filesystem
  phmap)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
