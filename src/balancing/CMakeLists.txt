# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Eigen3 REQUIRED QUIET)

add_library(hictk_balancing INTERFACE)
add_library(hictk::balancing ALIAS hictk_balancing)
target_sources(hictk_balancing INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/vc_impl.hpp)

target_include_directories(hictk_balancing INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(
        hictk_balancing
        INTERFACE hictk_project_warnings
        hictk_project_options
        hictk::common
        hictk::cooler
        hictk::pixel)

target_link_system_libraries(hictk_balancing INTERFACE Eigen3::Eigen)

install(TARGETS hictk_balancing INCLUDES)
