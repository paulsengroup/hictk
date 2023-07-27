# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(genomic_interval INTERFACE)
add_library(hictk::genomic_interval ALIAS genomic_interval)

target_sources(
  genomic_interval
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(genomic_interval INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                      "$<INSTALL_INTERFACE:include>")
target_link_libraries(genomic_interval INTERFACE hictk::chromosome hictk::common hictk::numeric)

include(GNUInstallDirs)
install(DIRECTORY include/hictk DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
