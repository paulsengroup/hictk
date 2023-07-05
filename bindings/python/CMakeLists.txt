# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.27)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

project(
        hictkpy
        LANGUAGES CXX
        HOMEPAGE_URL https://github.com/paulsengroup/hictk
        DESCRIPTION "High-performance C++ library to read and write .cool files.")

add_subdirectory(src)
