# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(tsl-hopscotch-map CONFIG REQUIRED)

add_library(hictk_variant INTERFACE)
add_library(hictk::variant ALIAS hictk_variant)

target_sources(hictk_variant INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/variant_buff_impl.hpp)
target_include_directories(hictk_variant INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(hictk_variant INTERFACE hictk_project_warnings hictk_project_options)

install(TARGETS hictk_variant INCLUDES)
