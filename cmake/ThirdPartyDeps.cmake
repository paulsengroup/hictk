# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)
FetchContent_Declare(
  _hictk_mimalloc
  URL
    "${CMAKE_CURRENT_SOURCE_DIR}/external/mimalloc-v3.1.5.tar.xz"
  URL_HASH SHA256=239b9196d6000942f351077f42be6170e5e4bd5dc7ee2642370d8222e4e0c4bb
  SYSTEM
)
if(HICTK_WITH_MIMALLOC)
  message(STATUS "Building with mimalloc support")
  set(MI_OVERRIDE ON)
  set(MI_BUILD_OBJECT ON)
  set(MI_BUILD_STATIC OFF)
  set(MI_BUILD_SHARED OFF)
  set(MI_BUILD_TESTS OFF)
  set(MI_TRACK_ASAN ${OPT_ENABLE_SANITIZER_ADDRESS})

  FetchContent_MakeAvailable(_hictk_mimalloc)
  # TODO this only works on unix-like systems
  add_library(hictk::mimalloc ALIAS mimalloc-obj)

  add_library(_hictk_mimalloc_override STATIC)
  add_library(hictk::mimalloc_override ALIAS _hictk_mimalloc_override)

  set(HICTK_MIMALLOC_OVERRIDE_CPP_FILE "${CMAKE_CURRENT_BINARY_DIR}/generated/mimalloc_override/mimalloc_override.cpp")
  file(WRITE "${HICTK_MIMALLOC_OVERRIDE_CPP_FILE}" "#include \"mimalloc-new-delete.h\"")
  target_sources(_hictk_mimalloc_override PRIVATE "${HICTK_MIMALLOC_OVERRIDE_CPP_FILE}")

  target_link_libraries(_hictk_mimalloc_override PRIVATE hictk::mimalloc)
  set_target_properties(
    _hictk_mimalloc_override
    PROPERTIES
      INTERFACE_LINK_LIBRARIES_DIRECT
        $<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:hictk::mimalloc>
  )
endif()
