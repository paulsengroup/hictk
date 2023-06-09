# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)

add_library(hic INTERFACE)
add_library(hictk::hic ALIAS hic)

target_sources(
  hic
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/block_cache_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/block_reader_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/filestream_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_reader_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/footer_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/header_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/index_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp)

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
  _hictk_fmt_tgt
  $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,libdeflate::libdeflate_shared,libdeflate::libdeflate_static>
  _hictk_phmap_tgt
  std::filesystem)
