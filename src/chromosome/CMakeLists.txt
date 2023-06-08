# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_chromosome INTERFACE)
add_library(hictk::chromosome ALIAS hictk_chromosome)

target_sources(hictk_chromosome INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/chromosome_impl.hpp)
target_include_directories(hictk_chromosome INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(hictk_chromosome INTERFACE hictk_project_warnings hictk_project_options hictk::common)

install(TARGETS hictk_chromosome INCLUDES)
