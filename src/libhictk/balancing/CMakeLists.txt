# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(bshoshany-thread-pool REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(span-lite REQUIRED QUIET)
find_package(zstd REQUIRED QUIET)

add_library(hictk_internal_balancing INTERFACE)
add_library(hictk::balancing ALIAS hictk_internal_balancing)

target_sources(
  hictk_internal_balancing
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/common.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/ice.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/methods.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/scale.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/sparse_matrix.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/vc.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/weights.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/impl/ice_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/impl/scale_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/impl/sparse_matrix_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/impl/vc_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/balancing/impl/weights_impl.hpp"
)

target_include_directories(
  hictk_internal_balancing
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_balancing
  INTERFACE
    hictk::common
    hictk::default_delete_zstd
    hictk::filestream
    hictk::pixel
)

target_link_system_libraries(
  hictk_internal_balancing
  INTERFACE
    bshoshany-thread-pool::bshoshany-thread-pool
    nonstd::span-lite
    phmap
    "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)

target_compile_definitions(hictk_internal_balancing INTERFACE span_FEATURE_MAKE_SPAN=1)
