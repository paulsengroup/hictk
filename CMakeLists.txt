# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...4.0)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# Not ideal to use this global variable, but necessary to make sure that tooling and projects use the same version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# strongly encouraged to enable this globally to avoid conflicts between -Wpedantic being enabled and -std=c++20 and
# -std=gnu++20 for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

set(ENABLE_DEVELOPER_MODE OFF CACHE BOOL "Enable 'developer mode'")

if(NOT CMAKE_FIND_PACKAGE_PREFER_CONFIG)
  set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ON)
endif()

set(HICTK_PROJECT_VERSION_MAJOR 2)
set(HICTK_PROJECT_VERSION_MINOR 1)
set(HICTK_PROJECT_VERSION_PATCH 4)
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

include(FetchContent)
FetchContent_Declare(
  _hictk_project_options
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.36.6.tar.xz"
  URL_HASH SHA256=b4df2a4107847248d4711ee48477cdf93c1a1c7d51a931937759eb9ffba67af8
  SYSTEM
)
FetchContent_MakeAvailable(_hictk_project_options)

get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their associated conan packages available. You can reduce
    # this list to only the configuration types you use, but only if one is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(
      CMAKE_CONFIGURATION_TYPES
      Debug
      Release
      RelWithDebInfo
      CACHE STRING
      "Enabled build types"
      FORCE
    )
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

include("${_hictk_project_options_SOURCE_DIR}/src/DynamicProjectOptions.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

# dynamic_project_options sets recommended defaults and provides user and developer modes and full GUI support for
# choosing options at configure time

# for more flexibility, look into project_options() macro

# Any default can be overridden set(<feature_name>_DEFAULT <value>) - set default for both user and developer modes
# set(<feature_name>_DEVELOPER_DEFAULT <value>) - set default for developer mode set(<feature_name>_USER_DEFAULT
# <value>) - set default for user mode

# Initialize project_options variable related to this project This overwrites `project_options` and sets
# `project_warnings` uncomment the options to enable them:

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
set(ENABLE_SANITIZER_LEAK_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_POINTER_COMPARE_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_COMPARE_DEVELOPER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_POINTER_SUBTRACT_DEVELOPER_DEFAULT ON)
set(ENABLE_SANITIZER_UNDEFINED_USER_DEFAULT OFF)
set(ENABLE_SANITIZER_UNDEFINED_DEVELOPER_DEFAULT ON)

dynamic_project_options(
  PREFIX
  hictk
  # CLANG_WARNINGS # Override the defaults for the CLANG warnings GCC_WARNINGS   # Override the defaults for the GCC
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

target_compile_features(hictk_project_options INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")

# Tweak fmt
target_compile_definitions(
  hictk_project_options
  INTERFACE
    FMT_HEADER_ONLY
    FMT_ENFORCE_COMPILE_STRING
    FMT_USE_FULL_CACHE_DRAGONBOX
)
# Tweak spdlog
target_compile_definitions(hictk_project_options INTERFACE SPDLOG_FMT_EXTERNAL)

if(WIN32)
  target_compile_definitions(
    hictk_project_options
    INTERFACE
      NOMINMAX
      _CRT_SECURE_NO_WARNINGS
  )
endif()

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
option(HICTK_BUILD_TOOLS "Build cli tools" ON)

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
