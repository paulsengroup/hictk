# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(bshoshany-thread-pool REQUIRED)
find_package(phmap REQUIRED)
find_package(span-lite REQUIRED)
find_package(zstd REQUIRED)

add_library(balancing INTERFACE)
add_library(hictk::balancing ALIAS balancing)

target_sources(
  balancing
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  balancing
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  balancing
  INTERFACE
    hictk::common
    hictk::filestream
    hictk::pixel
)

target_link_system_libraries(
  balancing
  INTERFACE
  bshoshany-thread-pool::bshoshany-thread-pool
  nonstd::span-lite
  phmap
  "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)

target_compile_definitions(balancing INTERFACE span_FEATURE_MAKE_SPAN=1)
