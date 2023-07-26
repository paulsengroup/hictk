# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

FetchContent_Declare(
  _hictk_cli11
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cli11-v2.3.2.tar.xz
  URL_HASH SHA256=009b7e7a29a4c1768760df470f288e79d746532e5f666776edafb52f18960685
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_fast_float
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fast_float-v5.2.0.tar.xz
  URL_HASH SHA256=4c46c081d2098d1d39f70a003e0ada92959b305c121addab60a92de1cfffaae2
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_fmt
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/fmt-v10.0.0.tar.xz
  URL_HASH SHA256=8570604ab8bc1c4cf70c3eecd278c88be3acf941373374c4908ddf9e7ae84288
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_highfive
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/highfive-v2.7.1.tar.xz
  URL_HASH SHA256=951596d3e85bbc8c6ea00cd73ee76e2af203dd29febdce827016378d2f0925e8
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_libdeflate
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/libdeflate-v1.18.tar.xz
  URL_HASH SHA256=f1e1e2432f9329a5f53939527afb46c417c843520bd526be7f777ab270eb65a0
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_phmap
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-hashmap-v1.3.11.tar.xz
  URL_HASH SHA256=f8f672e9fefdaa5fba555a77ff1037d9003401344dd651e71c98212e3eaad8cc
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_readerwriterqueue
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/readerwriterqueue-v1.0.6.tar.xz
  URL_HASH SHA256=332dc71267b625e0402515417f0fb63977354d233fc4b04b1f0ad319ad43110c
  EXCLUDE_FROM_ALL SYSTEM)

FetchContent_Declare(
  _hictk_spdlog
  URL ${CMAKE_CURRENT_SOURCE_DIR}/external/spdlog-v1.12.0.tar.xz
  URL_HASH SHA256=a37bb250032861c468716861f76aa97192a31107308aea8bf21cb0ad23e8693a
  EXCLUDE_FROM_ALL SYSTEM)

set(LIBDEFLATE_BUILD_SHARED_LIB ${BUILD_SHARED_LIBS})
set(LIBDEFLATE_BUILD_STATIC_LIB NOT ${BUILD_SHARED_LIBS})
set(LIBDEFLATE_COMPRESSION_SUPPORT OFF)
set(LIBDEFLATE_BUILD_GZIP OFF)
FetchContent_MakeAvailable(
  _hictk_libdeflate
  _hictk_phmap
  _hictk_project_options)

if(HICTK_BUILD_TOOLS)
  FetchContent_MakeAvailable(_hictk_cli11)
endif()

# Setup fast_float
FetchContent_GetProperties(_hictk_fast_float)
if(NOT _hictk_fast_float_POPULATED)
  FetchContent_Populate(_hictk_fast_float)
endif()

add_library(_hictk_fast_float_tgt INTERFACE)
target_include_directories(_hictk_fast_float_tgt SYSTEM INTERFACE ${_hictk_fast_float_SOURCE_DIR}/include)

# Setup fmt
FetchContent_GetProperties(_hictk_fmt)
if(NOT _hictk_fmt_POPULATED)
  FetchContent_Populate(_hictk_fmt)
endif()

add_library(_hictk_fmt_tgt INTERFACE)
target_include_directories(_hictk_fmt_tgt SYSTEM INTERFACE ${_hictk_fmt_SOURCE_DIR}/include)
target_compile_definitions(_hictk_fmt_tgt INTERFACE FMT_HEADER_ONLY FMT_ENFORCE_COMPILE_STRING)

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

add_library(_hictk_highfive_tgt INTERFACE)
target_include_directories(_hictk_highfive_tgt SYSTEM INTERFACE ${_hictk_highfive_SOURCE_DIR}/include)

# Setup parallel_hashmap
add_library(_hictk_phmap_tgt INTERFACE)
target_include_directories(_hictk_phmap_tgt SYSTEM INTERFACE ${_hictk_phmap_SOURCE_DIR})

# Setup fast_float
FetchContent_GetProperties(_hictk_readerwriterqueue)
if(NOT _hictk_readerwriterqueue_POPULATED)
  FetchContent_Populate(_hictk_readerwriterqueue)
endif()

add_library(_hictk_readerwriterqueue_tgt INTERFACE)
target_include_directories(_hictk_readerwriterqueue_tgt SYSTEM INTERFACE ${_hictk_readerwriterqueue_SOURCE_DIR})

# Setup spdlog
FetchContent_GetProperties(_hictk_spdlog)
if(NOT _hictk_spdlog_POPULATED)
  FetchContent_Populate(_hictk_spdlog)
endif()

add_library(_hictk_spdlog_tgt INTERFACE)
target_include_directories(_hictk_spdlog_tgt SYSTEM INTERFACE ${_hictk_spdlog_SOURCE_DIR}/include)
target_compile_definitions(_hictk_spdlog_tgt INTERFACE SPDLOG_FMT_EXTERNAL)

# Setup project_options
include(${_hictk_project_options_SOURCE_DIR}/Index.cmake)

# Disable clang-tidy for external projects
target_disable_clang_tidy(_hictk_highfive_tgt)
target_disable_clang_tidy(_hictk_fmt_tgt)
target_disable_clang_tidy(_hictk_phmap_tgt)
target_disable_clang_tidy(_hictk_spdlog_tgt)

if(LIBDEFLATE_BUILD_SHARED_LIB)
  target_disable_clang_tidy(libdeflate_shared)
endif()
if(LIBDEFLATE_BUILD_STATIC_LIB)
  target_disable_clang_tidy(libdeflate_static)
endif()
