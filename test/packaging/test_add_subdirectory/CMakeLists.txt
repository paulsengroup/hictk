# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.30)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

project(hictk_test_add_subdirectory)

set(HICTK_ENABLE_TESTING OFF)
set(HICTK_BUILD_EXAMPLES OFF)
set(HICTK_BUILD_BENCHMARKS OFF)
set(HICTK_BUILD_PYTHON_BINDINGS OFF)
set(HICTK_BUILD_TOOLS OFF)
add_subdirectory(hictk_root)

add_executable(hictk_test_add_subdirectory main.cpp)
target_link_libraries(hictk_test_add_subdirectory PRIVATE hictk::libhictk)
