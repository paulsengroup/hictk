# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.30)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

project(hictk_test_find_package)

find_package(hictk REQUIRED)

add_executable(hictk_test_find_package main.cpp)
target_link_libraries(hictk_test_find_package PRIVATE hictk::libhictk)
