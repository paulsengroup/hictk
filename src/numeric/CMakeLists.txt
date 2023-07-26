# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_numeric INTERFACE)
add_library(hictk::numeric ALIAS hictk_numeric)
target_sources(hictk_numeric INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/numeric_utils_impl.hpp)

target_include_directories(hictk_numeric INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(hictk_numeric INTERFACE hictk_project_warnings hictk_project_options)

target_link_system_libraries(hictk_numeric INTERFACE _hictk_fast_float)

install(TARGETS hictk_numeric INCLUDES)
