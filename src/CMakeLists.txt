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
add_subdirectory(variant)

if(HICTK_BUILD_TOOLS)
  message(STATUS "Building cli tools.")
  add_subdirectory(hictk)
endif()
