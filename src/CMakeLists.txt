# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(libhictk)

if(HICTK_BUILD_TOOLS)
  message(STATUS "Building cli tools.")
  add_subdirectory(hictk)
endif()
