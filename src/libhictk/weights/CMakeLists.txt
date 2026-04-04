# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap REQUIRED QUIET)

add_library(hictk_internal_weights STATIC)
add_library(hictk::weights ALIAS hictk_internal_weights)

target_sources(
  hictk_internal_weights
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/weights.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/weights_impl.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/weights.cpp"
)

target_include_directories(
  hictk_internal_weights
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_weights
  PUBLIC
    hictk::common
    hictk::pixel
)

target_link_system_libraries(
  hictk_internal_weights
  PUBLIC
    phmap
)
