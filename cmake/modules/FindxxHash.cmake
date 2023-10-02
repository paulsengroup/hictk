# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(PkgConfig)
pkg_check_modules(PC_xxHash QUIET xxHash)

find_path(
  xxHash_INCLUDE_DIR
  NAMES xxhash.h
  PATHS ${PC_xxHash_INCLUDE_DIRS}
  PATH_SUFFIXES include)
find_library(
  xxHash_LIBRARY
  NAMES xxhash
  PATHS ${PC_xxHash_LIBRARY_DIRS}
  PATH_SUFFIXES lib)

set(xxHash_VERSION ${PC_xxHash_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  xxHash
  FOUND_VAR xxHash_FOUND
  REQUIRED_VARS xxHash_LIBRARY xxHash_INCLUDE_DIR
  VERSION_VAR xxHash_VERSION)

if(xxHash_FOUND)
  set(xxHash_LIBRARIES ${xxHash_LIBRARY})
  set(xxHash_INCLUDE_DIRS ${xxHash_INCLUDE_DIR})
  set(xxHash_DEFINITIONS ${PC_xxHash_CFLAGS_OTHER})
endif()

if(xxHash_FOUND AND NOT TARGET xxHash::xxhash)
  add_library(xxHash::xxhash UNKNOWN IMPORTED)
  set_target_properties(
    xxHash::xxhash
    PROPERTIES IMPORTED_LOCATION "${xxHash_LIBRARY}"
               INTERFACE_COMPILE_OPTIONS "${PC_xxHash_CFLAGS_OTHER}"
               INTERFACE_INCLUDE_DIRECTORIES "${xxHash_INCLUDE_DIR}")
endif()

mark_as_advanced(xxHash_INCLUDE_DIR xxHash_LIBRARY)

set(xxHash_VERSION_STRING ${xxHash_VERSION})
