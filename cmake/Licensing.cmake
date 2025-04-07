# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

function(ConfigureLicensing license_file input_config_file output_config_file)
  set(PRE_CONFIGURE_FILE "${input_config_file}")
  set(POST_CONFIGURE_FILE "${output_config_file}")

  file(READ "${license_file}" HICTK_LICENSE)

  if(NOT WIN32)
    file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
  endif()

  configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
endfunction()
