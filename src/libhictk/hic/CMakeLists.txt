# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(FMT CONFIG REQUIRED)
find_package(libdeflate CONFIG REQUIRED)
find_package(phmap CONFIG REQUIRED)

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 QUIET)
  if(NOT Eigen3_FOUND)
    find_package(
      Eigen3
      CONFIG
      QUIET
      REQUIRED)
  endif()
endif()

add_library(hictk_hic INTERFACE)
add_library(hictk::hic ALIAS hictk_hic)

target_sources(
  hictk_hic
  PUBLIC FILE_SET
         HEADERS
         BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/include
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/block_cache_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/block_reader_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/file_reader_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/filestream_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/footer_cache_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/footer_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/header_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/hic_file_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/index_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/pixel_selector_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/weight_cache_impl.hpp)

target_include_directories(hictk_hic INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(
  hictk_hic
  INTERFACE hictk_project_options
            hictk_project_warnings
            hictk::balancing
            hictk::common
            hictk::chromosome
            hictk::reference)

target_link_system_libraries(
  hictk_hic
  INTERFACE
  $<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>
  fmt::fmt
  libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
  phmap
  std::filesystem)
