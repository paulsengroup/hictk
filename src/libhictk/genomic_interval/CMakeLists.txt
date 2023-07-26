# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_genomic_interval INTERFACE)
add_library(hictk::genomic_interval ALIAS hictk_genomic_interval)

target_sources(
  hictk_genomic_interval
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/genomic_interval_impl.hpp)
target_include_directories(hictk_genomic_interval INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(
  hictk_genomic_interval
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::chromosome
            hictk::common
            hictk::numeric)

install(TARGETS hictk_genomic_interval FILE_SET HEADERS EXCLUDE_FROM_ALL)
