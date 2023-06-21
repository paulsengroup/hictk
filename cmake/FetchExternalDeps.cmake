# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

FetchContent_Declare(
  _hictk_catch2
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/catch2-v3.3.2.tar.xz
  URL_HASH SHA256=05ac83d7c65d6ee4a499b8c98b32ef2cb30280d7691037981feffe498992804f
  SYSTEM)

FetchContent_Declare(
  _hictk_cli11
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cli11-v2.3.2.tar.xz
  URL_HASH SHA256=009b7e7a29a4c1768760df470f288e79d746532e5f666776edafb52f18960685
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  CLI11)

FetchContent_Declare(
  _hictk_fast_float
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fast_float-v5.2.0.tar.xz
  URL_HASH SHA256=4c46c081d2098d1d39f70a003e0ada92959b305c121addab60a92de1cfffaae2
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  FastFloat)

FetchContent_Declare(
  _hictk_fmt
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fmt-v10.0.0.tar.xz
  URL_HASH SHA256=8570604ab8bc1c4cf70c3eecd278c88be3acf941373374c4908ddf9e7ae84288
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  FMT)

FetchContent_Declare(
  _hictk_phmap
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-hashmap-v1.3.11.tar.xz
  URL_HASH SHA256=f8f672e9fefdaa5fba555a77ff1037d9003401344dd651e71c98212e3eaad8cc
  SYSTEM)

FetchContent_Declare(
  _hictk_project_options
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.29.0.tar.xz
  URL_HASH SHA256=ee2836af616d42e22c61048f4aedafd104fbfd97cf52bdd122a212d98777304c
  SYSTEM)

FetchContent_Declare(
  _hictk_readerwriterqueue
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/readerwriterqueue-v1.0.6.tar.xz
  URL_HASH SHA256=332dc71267b625e0402515417f0fb63977354d233fc4b04b1f0ad319ad43110c
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  readerwriterqueue)

FetchContent_Declare(
  _hictk_spanlite
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/spanlite-v0.10.3.tar.xz
  URL_HASH SHA256=3bfccb1a2f246da92303185df67de3e57de454e4aa5d861d3bfc81bff9771559
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  span-lite)

FetchContent_Declare(
  _hictk_spdlog
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog-v1.11.0.tar.xz
  URL_HASH SHA256=7bb89d5baba54638a2107291c40f2972428ac32a3c65609b2ffedb2d295ca1ad
  SYSTEM
  FIND_PACKAGE_ARGS
  NAMES
  spdlog)

set(HIGHFIVE_PARALLEL_HDF5 OFF)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN OFF)
set(HIGHFIVE_USE_XTENSOR OFF)
set(HIGHFIVE_USE_OPENCV OFF)

FetchContent_MakeAvailable(
  _hictk_catch2
  _hictk_cli11
  _hictk_fast_float
  _hictk_fmt
  _hictk_phmap
  _hictk_project_options
  _hictk_readerwriterqueue
  _hictk_spanlite
  _hictk_spdlog)

include_directories(SYSTEM ${_hictk_phmap_SOURCE_DIR})
include(${_hictk_project_options_SOURCE_DIR}/Index.cmake)

list(APPEND CMAKE_MODULE_PATH ${_hictk_catch2_SOURCE_DIR}/extras)
