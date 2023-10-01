# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(HICTK_EXAMPLES quicktour_example_001 quicktour_example_002 quicktour_example_003)

foreach(EXAMPLE ${HICTK_EXAMPLES})
  add_executable("hictk_${EXAMPLE}" "${EXAMPLE}.cpp")
  target_link_libraries(
    "hictk_${EXAMPLE}"
    PRIVATE hictk_project_warnings
            hictk_project_options
            hictk::file
            hictk::transformers)
endforeach()
