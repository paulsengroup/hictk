# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(HICTK_WITH_ARROW)
  find_package(Arrow REQUIRED QUIET)
endif()

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 REQUIRED QUIET)
endif()

find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_transformers INTERFACE)
add_library(hictk::transformers ALIAS hictk_internal_transformers)

target_sources(
  hictk_internal_transformers
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
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
    hictk::cooler
    hictk::file
    hictk::hic
)

target_link_system_libraries(hictk_internal_transformers INTERFACE "$<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>" spdlog::spdlog_header_only)

if(HICTK_WITH_ARROW_SHARED)
  target_link_system_libraries(hictk_internal_transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_shared>")
else()
  target_link_system_libraries(hictk_internal_transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_static>")
endif()
