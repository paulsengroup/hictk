# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(bshoshany-thread-pool REQUIRED QUIET)
find_package(concurrentqueue REQUIRED QUIET)
find_package(fmt REQUIRED QUIET)
find_package(libdeflate REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(readerwriterqueue REQUIRED QUIET)
find_package(zstd REQUIRED QUIET)

set(HICTK_CXX_BYTE_ORDER "${CMAKE_CXX_BYTE_ORDER}" CACHE STRING "Specify the CPU architecture endianness")

if(HICTK_CXX_BYTE_ORDER STREQUAL "BIG_ENDIAN")
  message(STATUS "Detected ${HICTK_CXX_BYTE_ORDER} architecture")
  message(FATAL_ERROR "hictk::hic cannot be built on machines with big-endian architecture.")
elseif(HICTK_CXX_BYTE_ORDER STREQUAL "LITTLE_ENDIAN")
  message(STATUS "Detected ${HICTK_CXX_BYTE_ORDER} architecture")
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
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/block_reader.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/cache.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/common.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/file_reader.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/file_writer.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/file_writer_data_structures.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/file_zoomify.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/footer.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/footer_cache.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/header.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/index.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/interaction_block.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/interaction_to_block_mapper.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/pixel_selector.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/serialized_block_pqueue.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/utils.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/validation.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/block_cache_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/block_reader_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/file_reader_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/file_writer_data_structures_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/file_writer_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/file_zoomify_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/footer_cache_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/footer_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/header_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/hic_file_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/index_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/interaction_block_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/interaction_to_block_mapper_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/pixel_selector_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/serialized_block_pqueue_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/utils_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/utils_merge_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/validation_impl.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/hic/impl/weight_cache_impl.hpp"
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
    hictk::default_delete_libdeflate
    hictk::default_delete_zstd
    hictk::binary_buffer
    hictk::chromosome
    hictk::expected_values_aggregator
    hictk::filestream
    hictk::hash
    hictk::pixel
    hictk::reference
    hictk::static_binary_buffer
    hictk::tmpdir
    hictk::transformers
    hictk::version
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
