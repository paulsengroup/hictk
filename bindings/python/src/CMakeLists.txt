# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(
  Python 3.7
  COMPONENTS Interpreter Development
  REQUIRED)
find_package(pybind11 REQUIRED CONFIG)

#find_library(ARROW_LIBS arrow HINTS ${HICTK_ARROW_ROOT} NO_DEFAULT_PATH)
find_library(PYARROW_LIBS arrow_python HINTS ${HICTK_ARROW_ROOT} NO_DEFAULT_PATH)

#if (ARROW_LIBS STREQUAL "ARROW_LIBS-NOTFOUND")
#  message(FATAL_ERROR "Unable to find libarrow under ${HICTK_ARROW_ROOT}")
#endif()
if(PYARROW_LIBS STREQUAL "PYARROW_LIBS-NOTFOUND")
  message(FATAL_ERROR "Unable to find libarrow_python under ${HICTK_ARROW_ROOT}")
endif()

pybind11_add_module(hictkpy MODULE hictk.cpp)

target_link_directories(hictkpy PUBLIC "${HICTK_ARROW_ROOT}")
target_include_directories(hictkpy PRIVATE "${HICTK_ARROW_ROOT}/include")
target_link_libraries(
  hictkpy
  PRIVATE hictk_project_options
          hictk_project_warnings
  PUBLIC  hictk::cooler
          hictk::hic)

target_link_system_libraries(
        hictkpy
        PUBLIC
        arrow_python)
        #${ARROW_LIBS})
