# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(HICTK_LICENSE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(PRE_CONFIGURE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/libhictk/config/license.hpp.in")
set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/libhictk/common/version/include/hictk/license.hpp")
file(READ "${HICTK_LICENSE_FILE}" HICTK_LICENSE)

if(NOT WIN32)
  file(LOCK "${POST_CONFIGURE_FILE}" GUARD FILE)
endif()

configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
