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

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 QUIET REQUIRED)
endif()

add_library(hic INTERFACE)
add_library(hictk::hic ALIAS hic)

target_sources(
  hic
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(hic INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                         "$<INSTALL_INTERFACE:include>")
target_link_libraries(
  hic
  INTERFACE hictk::balancing
            hictk::common
            hictk::binary_buffer
            hictk::chromosome
            hictk::expected_values_aggregator
            hictk::filestream
            hictk::pixel
            hictk::reference
            hictk::transformers)

target_link_system_libraries(
  hic
  INTERFACE
  bshoshany-thread-pool::bshoshany-thread-pool
  concurrentqueue::concurrentqueue
  "$<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>"
  fmt::fmt-header-only
  "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
  phmap
  readerwriterqueue::readerwriterqueue
  "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>")
