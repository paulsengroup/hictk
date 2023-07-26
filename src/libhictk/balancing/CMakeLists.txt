# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap CONFIG REQUIRED)

add_library(balancing INTERFACE)
add_library(hictk::balancing ALIAS balancing)

target_sources(
  balancing
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(balancing INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                               $<INSTALL_INTERFACE:include>)
target_link_libraries(balancing INTERFACE hictk::common hictk::pixel)

target_link_system_libraries(balancing INTERFACE phmap)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
