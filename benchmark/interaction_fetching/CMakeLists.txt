# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

target_sources(
  hictk_benchmark
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler_gw_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cooler_trans_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_gw_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/file_trans_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_cis_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_gw_queries.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic_trans_queries.cpp"
)

target_link_libraries(
  hictk_benchmark
  PUBLIC
    hictk::cooler
    hictk::file
    hictk::hic
)
