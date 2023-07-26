# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_transformers INTERFACE)
add_library(hictk::transformers ALIAS hictk_transformers)

target_sources(
  hictk_transformers
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(hictk_transformers INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                                        $<INSTALL_INTERFACE:include>)
target_link_libraries(hictk_transformers INTERFACE hictk::cooler hictk::hic)

include(GNUInstallDirs)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
