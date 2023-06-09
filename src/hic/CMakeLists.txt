# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt QUIET REQUIRED)
find_package(Filesystem REQUIRED)
find_package(libdeflate QUIET REQUIRED)
find_package(tsl-ordered-map QUIET REQUIRED)

add_library(hic INTERFACE)
add_library(hictk::hic ALIAS hic)

target_sources(
  hic
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cache_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/filestream_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_stream_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_utils_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_footer_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_header_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_matrix_selector_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/filestream_impl.hpp)

target_include_directories(hic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
  hic
  INTERFACE hictk_project_options
            hictk_project_warnings
            hictk::common
            hictk::chromosome
            hictk::reference)

target_link_system_libraries(
  hic
  INTERFACE
  fmt::fmt
  libdeflate::libdeflate
  tsl::ordered_map
  std::filesystem)
