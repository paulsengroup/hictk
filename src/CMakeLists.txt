# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(libhictk)

if(HICTK_BUILD_TOOLS)
  add_subdirectory(hictk)
endif()
