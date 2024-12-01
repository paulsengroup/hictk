# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

target_sources(
  hictk_benchmark
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/bin_ids_to_genomic_coords.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/genomic_coords_to_bin_ids.cpp"
)

target_link_libraries(hictk_benchmark PUBLIC hictk::bin_table)
