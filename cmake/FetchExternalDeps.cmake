# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

FetchContent_MakeAvailable(_hictk_project_options)

include(${_hictk_project_options_SOURCE_DIR}/Index.cmake)
