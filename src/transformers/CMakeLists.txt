# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_transformers INTERFACE)
add_library(hictk::transformers ALIAS hictk_transformers)

target_sources(hictk_transformers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/transformers_impl.hpp)
target_include_directories(hictk_transformers INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(
        hictk_transformers
        INTERFACE hictk_project_warnings
        hictk_project_options
        hictk::cooler
        hictk::hic)

install(TARGETS hictk_transformers INCLUDES)

