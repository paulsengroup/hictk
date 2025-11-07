# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...4.1)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

set(ENABLE_DEVELOPER_MODE OFF CACHE BOOL "Enable 'developer mode'")

if(NOT DEFINED CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

set(HICTK_PROJECT_VERSION_MAJOR 2)
set(HICTK_PROJECT_VERSION_MINOR 1)
set(HICTK_PROJECT_VERSION_PATCH 5)
set(HICTK_PROJECT_VERSION_SUFFIX "")

project(
  hictk
  LANGUAGES
    C
    CXX
  VERSION "${HICTK_PROJECT_VERSION_MAJOR}.${HICTK_PROJECT_VERSION_MINOR}.${HICTK_PROJECT_VERSION_PATCH}"
  HOMEPAGE_URL https://github.com/paulsengroup/hictk
  DESCRIPTION "Blazing fast toolkit to work with .hic and .cool files."
)

if(DEFINED CMAKE_CXX_STANDARD AND CMAKE_CXX_STANDARD LESS 17)
  message(
    FATAL_ERROR
    "hictk requires C++17 or newer (found CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}): please re-configure the project using -DCMAKE_CXX_STANDARD=17"
  )
endif()

include(FetchContent)
FetchContent_Declare(
  _hictk_project_options
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.36.6.tar.xz"
  URL_HASH SHA256=b4df2a4107847248d4711ee48477cdf93c1a1c7d51a931937759eb9ffba67af8
  SYSTEM
)
FetchContent_MakeAvailable(_hictk_project_options)

include("${_hictk_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

set(ENABLE_CACHE_DEFAULT ON)
set(ENABLE_COMPILE_COMMANDS_SYMLINK_DEFAULT OFF)
set(ENABLE_CONAN_DEFAULT OFF)
set(ENABLE_CPPCHECK_DEFAULT OFF)
set(ENABLE_DOXYGEN_DEFAULT OFF)
set(ENABLE_INTERPROCEDURAL_OPTIMIZATION_DEFAULT ON)
set(ENABLE_NATIVE_OPTIMIZATION_DEFAULT OFF)
set(ENABLE_PCH_DEFAULT OFF)

set(ENABLE_SANITIZER_ADDRESS_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_ADDRESS_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_LEAK_USER_DEFAULT OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
  set(ENABLE_SANITIZER_LEAK_DEVELOPER_DEFAULT ON)
else()
  set(ENABLE_SANITIZER_LEAK_DEVELOPER_DEFAULT OFF)
endif()
set(ENABLE_SANITIZER_POINTER_COMPARE_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_COMPARE_DEVELOPER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_UNDEFINED_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_UNDEFINED_DEVELOPER_DEFAULT ON)

dynamic_project_options(
  PREFIX
  hictk
  CPPCHECK_OPTIONS
  --enable=performance,portability,style,warning
  --inline-suppr
  # We cannot act on a bug/missing feature of cppcheck
  --suppress=internalAstError
  # if a file does not have an internalAstError, we get an unmatchedSuppression error
  --suppress=unmatchedSuppression
  --suppress=passedByValue
  --inconclusive
  MSVC_WARNINGS
  "${MSVC_WARNINGS}"
  CLANG_WARNINGS
  "${CLANG_WARNINGS}"
  GCC_WARNINGS
  "${GCC_WARNINGS}"
  CUDA_WARNINGS
  "${CUDA_WARNINGS}"
)

add_library(hictk::project_options ALIAS hictk_project_options)
add_library(hictk::project_warnings ALIAS hictk_project_warnings)

target_compile_features(
  hictk_project_options
  INTERFACE
    "cxx_std_$<IF:$<BOOL:${CMAKE_CXX_STANDARD}>,17,${CMAKE_CXX_STANDARD}>"
)

target_compile_definitions(
  hictk_project_options
  INTERFACE
    # Tweak fmt
    FMT_HEADER_ONLY
    # starting with fmt v12, clang-tidy breaks when FMT_ENFORCE_COMPILE_STRING is defined
    # https://github.com/fmtlib/fmt/commit/619b3a5aa031a25e8d06fec1e03992e5e727c5a5#diff-bdc6f79e8e9f5b4331d66fb785636a87d29f55cf729865e13925b4209424c878R4223
    $<$<NOT:$<BOOL:OPT_ENABLE_CLANG_TIDY>>:FMT_ENFORCE_COMPILE_STRING>
    FMT_USE_FULL_CACHE_DRAGONBOX
    # Tweak spdlog
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_$<IF:$<CONFIG:Debug>,DEBUG,INFO>
    SPDLOG_CLOCK_COARSE
    SPDLOG_FMT_EXTERNAL
    SPDLOG_NO_THREAD_ID
    SPDLOG_NO_ATOMIC_LEVELS
    # Windows-specific tweaks
    $<$<BOOL:${WIN32}>:NOMINMAX>
    $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>
)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif()

include(CMakeDependentOption)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
option(HICTK_ENABLE_TESTING "Build unit tests" ON)
option(HICTK_ENABLE_FUZZY_TESTING "Build fuzzy tests" OFF)
option(HICTK_ENABLE_GIT_VERSION_TRACKING "Retrieve project version and metadata from git" ON)
option(HICTK_BUILD_EXAMPLES "Build examples" OFF)
option(HICTK_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(HICTK_WITH_ARROW "Build with arrow support" ON)
option(HICTK_WITH_ARROW_SHARED "Force dynamic linking to Arrow libs" OFF)
option(HICTK_WITH_EIGEN "Build with Eigen3 support" ON)
option(HICTK_WITH_MIMALLOC "Statically link mimalloc to hictk to attempt to override malloc and new" OFF)
option(HICTK_BUILD_TOOLS "Build cli tools" ON)
option(
  HICTK_ASSUME_STD_FILESYSTEM_AVAIL
  "Assume that the <filesystem> header is available and does not require special linker flags"
  OFF
)

cmake_dependent_option(
  HICTK_DOWNLOAD_TEST_DATASET
  "Download datasets required by unit and integration tests"
  ON
  "HICTK_ENABLE_TESTING"
  OFF
)
cmake_dependent_option(
  HICTK_ENABLE_TELEMETRY
  "Build CLI tools with support for telemetry"
  ON
  "HICTK_BUILD_TOOLS"
  OFF
)

set(
  HICTK_EXPORTER_OTLP_ENDPOINT
  # "http://localhost:9"
  "https://hictk-telemetry.paulsenlab.com:4318"
  CACHE STRING
  "Endpoint where to submit telemetry data"
)

include(cmake/ThirdPartyDeps.cmake)

if(NOT HICTK_ENABLE_TELEMETRY)
  unset(HICTK_EXPORTER_OTLP_ENDPOINT CACHE)
endif()

if(HICTK_WITH_EIGEN)
  target_compile_definitions(hictk_project_options INTERFACE HICTK_WITH_EIGEN)
endif()

if(HICTK_WITH_ARROW)
  target_compile_definitions(hictk_project_options INTERFACE HICTK_WITH_ARROW)
endif()

if(BUILD_SHARED_LIBS)
  set(HICTK_WITH_ARROW_SHARED ON CACHE BOOL "Force dynamic linking to Arrow libs" FORCE)
endif()

if(HICTK_ASSUME_STD_FILESYSTEM_AVAIL)
  add_library(hictk_internal_std_filesystem_dummy INTERFACE)
  add_library(std::filesystem ALIAS hictk_internal_std_filesystem_dummy)
endif()

add_subdirectory(src)

if(HICTK_ENABLE_TESTING)
  enable_testing()
  target_compile_definitions(hictk_project_options INTERFACE HICTK_ENABLE_TESTING)
  add_subdirectory(test)
endif()

if(HICTK_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(HICTK_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()

include(cmake/Install.cmake)
