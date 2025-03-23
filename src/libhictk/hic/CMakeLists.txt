# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(bshoshany-thread-pool REQUIRED)
find_package(concurrentqueue REQUIRED)
find_package(FMT REQUIRED)
find_package(libdeflate REQUIRED)
find_package(phmap REQUIRED)
find_package(readerwriterqueue REQUIRED)
find_package(zstd REQUIRED)

set(HICTK_CXX_BYTE_ORDER "${CMAKE_CXX_BYTE_ORDER}" CACHE STRING "Specify the CPU architecture endianness")

if(HICTK_CXX_BYTE_ORDER STREQUAL "BIG_ENDIAN")
  message(STATUS "Detected ${HICTK_CXX_BYTE_ORDER} architecture.")
  message(FATAL_ERROR "hictk::hic cannot be built on machines with big-endian architecture.")
elseif(HICTK_CXX_BYTE_ORDER STREQUAL "LITTLE_ENDIAN")
  message(STATUS "Detected ${HICTK_CXX_BYTE_ORDER} architecture.")
else()
  message(
    FATAL_ERROR
    "Unable to detect the CPU architecture endianness.\nPlease specify the endianness of the target architecture with -DHICTK_CXX_BYTE_ORDER=*_ENDIAN.\nSee https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_BYTE_ORDER.html for more details."
  )
endif()

add_library(hictk_internal_hic INTERFACE)
add_library(hictk::hic ALIAS hictk_internal_hic)

target_sources(
  hictk_internal_hic
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_hic
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_hic
  INTERFACE
    hictk::balancing
    hictk::common
    hictk::binary_buffer
    hictk::chromosome
    hictk::expected_values_aggregator
    hictk::filestream
    hictk::pixel
    hictk::reference
    hictk::transformers
)

target_link_system_libraries(
  hictk_internal_hic
  INTERFACE
  bshoshany-thread-pool::bshoshany-thread-pool
  concurrentqueue::concurrentqueue
  fmt::fmt-header-only
  "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
  phmap
  readerwriterqueue::readerwriterqueue
  "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)
