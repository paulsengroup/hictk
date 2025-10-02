# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(HICTK_WITH_ARROW)
  find_package(Arrow REQUIRED QUIET)
endif()

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 REQUIRED QUIET)
endif()

find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_transformers INTERFACE)
add_library(hictk::transformers ALIAS hictk_internal_transformers)

target_sources(
  hictk_internal_transformers
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/coarsen.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/common.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/diagonal_band.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/join_genomic_coords.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/pixel_merger.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/stats.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/to_dataframe.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/to_dense_matrix.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/to_sparse_matrix.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/coarsen_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/diagonal_band_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/join_genomic_coords_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/pixel_merger_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/stats_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/to_dataframe_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/to_dense_matrix_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/transformers/impl/to_sparse_matrix_impl.hpp"
)

target_include_directories(
  hictk_internal_transformers
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_transformers
  INTERFACE
    hictk::common
    hictk::cooler
    hictk::file
    hictk::hic
    hictk::type_traits
)

target_link_system_libraries(
  hictk_internal_transformers
  INTERFACE
    "$<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>"
    spdlog::spdlog_header_only
    phmap
)

if(HICTK_WITH_ARROW_SHARED)
  target_link_system_libraries(hictk_internal_transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_shared>")
else()
  target_link_system_libraries(hictk_internal_transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_static>")
endif()
