# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)
find_package(libdeflate REQUIRED)
find_package(zstd REQUIRED)

add_library(common INTERFACE)
add_library(hictk::common ALIAS common)

target_sources(
  common
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(common INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                            "$<INSTALL_INTERFACE:include>")

target_link_system_libraries(common INTERFACE fmt::fmt-header-only)

if(BUILD_SHARED_LIBS AND HICTK_PREFER_STATIC_LINKING)
  if(TARGET libdeflate::libdeflate_static)
    target_link_system_libraries(common INTERFACE libdeflate::libdeflate_static)
  else()
    target_link_system_libraries(common INTERFACE libdeflate::libdeflate_shared)
  endif()
  if(TARGET zstd::libzstd_static)
    target_link_system_libraries(common INTERFACE zstd::libzstd_static)
  else()
    target_link_system_libraries(common INTERFACE zstd::libzstd_shared)
  endif()
else()
  target_link_system_libraries(
    common
    INTERFACE
    "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
    "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>")
endif()
