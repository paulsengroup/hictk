# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include("${PROJECT_SOURCE_DIR}/cmake/Licensing.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/Versioning.cmake")

set(HICTK_VERSION_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/include")

# gersemi: off
ConfigureLicensing(
  "${PROJECT_SOURCE_DIR}/LICENSE"
  "${CMAKE_CURRENT_SOURCE_DIR}/license.hpp.in"
  "${HICTK_VERSION_INCLUDE_DIR}/hictk/license.hpp"
)

ConfigureVersioning(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${HICTK_VERSION_INCLUDE_DIR}/hictk"
)
# gersemi: on

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_common_version INTERFACE)
add_library(hictk::version ALIAS hictk_internal_common_version)
add_dependencies(hictk_internal_common_version _hictk_check_git)

target_sources(
  hictk_internal_common_version
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${HICTK_VERSION_INCLUDE_DIR}"
      FILES
        "${HICTK_VERSION_INCLUDE_DIR}/hictk/git.hpp"
        "${HICTK_VERSION_INCLUDE_DIR}/hictk/license.hpp"
        "${HICTK_VERSION_INCLUDE_DIR}/hictk/version.hpp"
)

target_include_directories(
  hictk_internal_common_version
  INTERFACE
    "$<BUILD_INTERFACE:${HICTK_VERSION_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(hictk_internal_common_version INTERFACE fmt::fmt-header-only)
