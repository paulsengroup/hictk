# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.7
  COMPONENTS Interpreter Development
  REQUIRED)

pybind11_add_module(hictkpy MODULE hictkpy.cpp)

target_link_libraries(
  hictkpy
  PRIVATE hictk_project_options hictk_project_warnings
  PUBLIC hictk::cooler hictk::hic)
