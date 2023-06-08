# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(HICXX_USE_CURL)
    find_package(CURL REQUIRED)
endif()

find_package(libdeflate QUIET REQUIRED)

find_package(fmt QUIET REQUIRED)

find_package(Filesystem REQUIRED)

find_package(tsl-ordered-map QUIET REQUIRED)

add_library(hicxx INTERFACE)

target_sources(
    hicxx
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cache_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_utils_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_footer_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_header_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/hic_matrix_selector_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/local_fstream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/remote_fstream_impl.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/hicxx.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/cache.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/common.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_footer.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_matrix_selector.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_file_stream.hpp
              ${CMAKE_CURRENT_SOURCE_DIR}/include/hicxx/internal/hic_header.hpp)

target_include_directories(hicxx INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
    hicxx
    INTERFACE hicxx_project_options
              hicxx_project_warnings
              $<$<BOOL:${HICXX_USE_CURL}>:CURL::libcurl>
              fmt::fmt
              libdeflate::libdeflate
              tsl::ordered_map
              std::filesystem)

target_compile_definitions(hicxx INTERFACE FMT_ENFORCE_COMPILE_STRING $<$<BOOL:${HICXX_USE_CURL}>:HICXX_USE_CURL>)
