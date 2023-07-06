# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

FetchContent_Declare(
  _hictk_cli11
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cli11-v2.3.2.tar.xz
  URL_HASH SHA256=009b7e7a29a4c1768760df470f288e79d746532e5f666776edafb52f18960685
  SYSTEM)

FetchContent_Declare(
  _hictk_fast_float
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fast_float-v5.2.0.tar.xz
  URL_HASH SHA256=4c46c081d2098d1d39f70a003e0ada92959b305c121addab60a92de1cfffaae2
  SYSTEM)

FetchContent_Declare(
  _hictk_fmt
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fmt-v9.1.0.tar.xz
  URL_HASH SHA256=d2b242c76dbd3c7e0d763cb8b9021887e4b4e04f4adada24b0f19d4edbf02f96
  SYSTEM)

FetchContent_Declare(
  _hictk_highfive
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/highfive-v2.7.1.tar.xz
  URL_HASH SHA256=951596d3e85bbc8c6ea00cd73ee76e2af203dd29febdce827016378d2f0925e8
  SYSTEM)

FetchContent_Declare(
  _hictk_libdeflate
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libdeflate-v1.18.tar.xz
  URL_HASH SHA256=f1e1e2432f9329a5f53939527afb46c417c843520bd526be7f777ab270eb65a0
  FIND_PACKAGE_ARGS
  NAMES
  libdeflate
  SYSTEM)

FetchContent_Declare(
  _hictk_phmap
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-hashmap-v1.3.11.tar.xz
  URL_HASH SHA256=f8f672e9fefdaa5fba555a77ff1037d9003401344dd651e71c98212e3eaad8cc
  SYSTEM)

FetchContent_Declare(
  _hictk_readerwriterqueue
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/readerwriterqueue-v1.0.6.tar.xz
  URL_HASH SHA256=332dc71267b625e0402515417f0fb63977354d233fc4b04b1f0ad319ad43110c
  SYSTEM)

FetchContent_Declare(
  _hictk_spdlog
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog-v1.11.0.tar.xz
  URL_HASH SHA256=7bb89d5baba54638a2107291c40f2972428ac32a3c65609b2ffedb2d295ca1ad
  FIND_PACKAGE_ARGS
  NAMES
  spdlog
  VERSION
  1.11
  SYSTEM)

FetchContent_MakeAvailable(
  _hictk_fast_float
  _hictk_libdeflate
  _hictk_phmap
  _hictk_project_options)

if(HICTK_BUILD_TOOLS)
  FetchContent_MakeAvailable(
    _hictk_cli11
    _hictk_fast_float
    _hictk_readerwriterqueue
    _hictk_spdlog)
endif()

# Setup HighFive
set(HIGHFIVE_PARALLEL_HDF5 OFF)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN OFF)
set(HIGHFIVE_USE_XTENSOR OFF)
set(HIGHFIVE_USE_OPENCV OFF)
FetchContent_GetProperties(_hictk_highfive)
if(NOT _hictk_highfive_POPULATED)
  FetchContent_Populate(_hictk_highfive)
endif()

add_library(HighFive INTERFACE)
target_include_directories(HighFive INTERFACE ${_hictk_highfive_SOURCE_DIR}/include)

# Setup fmt
FetchContent_GetProperties(_hictk_fmt)
if(NOT _hictk_fmt_POPULATED)
  FetchContent_Populate(_hictk_fmt)
endif()

add_library(_hictk_fmt_tgt INTERFACE)
target_include_directories(_hictk_fmt_tgt INTERFACE ${_hictk_fmt_SOURCE_DIR}/include)

# Setup parallel_hashmap
add_library(_hictk_phmap_tgt INTERFACE)
target_include_directories(_hictk_phmap_tgt INTERFACE ${_hictk_phmap_SOURCE_DIR})

# Setup project_options
include(${_hictk_project_options_SOURCE_DIR}/Index.cmake)
