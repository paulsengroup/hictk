# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_common INTERFACE)
add_library(hictk::common ALIAS hictk_common)

target_sources(
  hictk_common
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(hictk_common INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                                  $<INSTALL_INTERFACE:include>)

include(GNUInstallDirs)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
