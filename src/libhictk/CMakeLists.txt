# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(balancing)
add_subdirectory(bin_table)
add_subdirectory(chromosome)
add_subdirectory(common)
add_subdirectory(cooler)
add_subdirectory(expected_values_aggregator)
add_subdirectory(file)
add_subdirectory(filestream)
add_subdirectory(formatting)
add_subdirectory(genomic_interval)
add_subdirectory(hic)
add_subdirectory(numeric)
add_subdirectory(pixel)
add_subdirectory(reference)
add_subdirectory(transformers)
add_subdirectory(variant)

add_library(libhictk INTERFACE)
add_library(hictk::libhictk ALIAS libhictk)

target_link_libraries(
  libhictk
  INTERFACE hictk::balancing
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::cooler
            hictk::expected_values_aggregator
            hictk::file
            hictk::filestream
            hictk::format
            hictk::genomic_interval
            hictk::hic
            hictk::numeric
            hictk::pixel
            hictk::reference
            hictk::transformers
            hictk::variant)

target_include_directories(
  libhictk
  INTERFACE
  INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" "$<INSTALL_INTERFACE:include>")
