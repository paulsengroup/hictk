# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(balancing)
add_subdirectory(bin_table)
add_subdirectory(chromosome)
add_subdirectory(common)
add_subdirectory(cooler)
add_subdirectory(fmt)
add_subdirectory(genomic_interval)
add_subdirectory(hic)
add_subdirectory(numeric)
add_subdirectory(pixel)
add_subdirectory(reference)
add_subdirectory(transformers)
add_subdirectory(variant)

add_library(libhictk INTERFACE)
add_library(hictk::hictk ALIAS libhictk)

target_link_libraries(
  libhictk
  INTERFACE hictk::balancing
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::cooler
            hictk::fmt
            hictk::genomic_interval
            hictk::hic
            hictk::numeric
            hictk::pixel
            hictk::reference
            hictk::transformers
            hictk::variant)
