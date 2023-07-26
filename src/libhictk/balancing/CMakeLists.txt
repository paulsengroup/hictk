# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap CONFIG REQUIRED)

add_library(hictk_balancing INTERFACE)
add_library(hictk::balancing ALIAS hictk_balancing)

target_sources(
  hictk_balancing
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(hictk_balancing INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                                     $<INSTALL_INTERFACE:include>)
target_link_libraries(hictk_balancing INTERFACE hictk::common hictk::pixel)

target_link_system_libraries(hictk_balancing INTERFACE phmap)

include(GNUInstallDirs)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
