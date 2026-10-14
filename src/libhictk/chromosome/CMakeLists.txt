# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_chromosome STATIC)
add_library(hictk::chromosome ALIAS hictk_internal_chromosome)

target_sources(
  hictk_internal_chromosome
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/chromosome.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/chromosome_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/chromosome.cpp"
)

target_include_directories(
  hictk_internal_chromosome
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(hictk_internal_chromosome PRIVATE hictk::hash)
