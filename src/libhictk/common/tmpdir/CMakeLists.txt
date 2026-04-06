# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_common_tmpdir STATIC)
add_library(hictk::tmpdir ALIAS hictk_internal_common_tmpdir)

include(CheckCXXSourceCompiles)

message(STATUS "Checking if mkdtemp() is available...")

check_cxx_source_compiles(
  [[
#include <unistd.h>

#include <cstdio>
#include <string>


int main() {
    std::string dir{"/tmp/hictk-tmp-XXXXXXXXXX"};
    if (!mkdtemp(dir.data())) {
        return 1;
    }
    return 0;
}
]]
  HICTK_MKDTEMP_AVAIL
)

if(HICTK_MKDTEMP_AVAIL)
  set(HICTK_MKDTEMP_AVAIL true)
else()
  set(HICTK_MKDTEMP_AVAIL false)
endif()

message(STATUS "mkdtemp() available: ${HICTK_MKDTEMP_AVAIL}")

target_sources(
  hictk_internal_common_tmpdir
  PUBLIC
    FILE_SET HEADERS
      BASE_DIRS
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/tmpdir.hpp"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/tmpdir.cpp"
)

target_compile_definitions(hictk_internal_common_tmpdir PRIVATE $<$<BOOL:${HICTK_MKDTEMP_AVAIL}>:HICTK_MKDTEMP_AVAIL>)

target_include_directories(
  hictk_internal_common_tmpdir
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_common_tmpdir
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
)
