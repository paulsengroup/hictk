# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

target_sources(hictk_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/fetch_chromosomes.cpp")

target_link_libraries(hictk_benchmark PUBLIC hictk::reference)
