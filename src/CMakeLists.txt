# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(libhictk)

if(HICTK_BUILD_TOOLS)
  option(HICTK_ENABLE_TELEMETRY "Build CLI tools with support for telemetry" ON)
  message(STATUS "Building CLI tools.")
  add_subdirectory(hictk)
endif()
