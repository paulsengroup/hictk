# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_tools_build_metadata STATIC)

target_sources(
  hictk_internal_tools_build_metadata
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/build_options.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/dependency_metadata.cpp"
)

target_include_directories(hictk_internal_tools_build_metadata PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_definitions(
  hictk_internal_tools_build_metadata
  PUBLIC
    HICTK_SYSTEM_NAME=\"${CMAKE_SYSTEM_NAME}\"
    HICTK_SYSTEM_VERSION=\"${CMAKE_HOST_SYSTEM_VERSION}\"
    HICTK_SYSTEM_PROCESSOR=\"${CMAKE_SYSTEM_PROCESSOR}\"
    HICTK_BUILD_TYPE=\"$<CONFIG>\"
    HICTK_CXX_COMPILER_ID=\"$<CXX_COMPILER_ID>\"
    HICTK_CXX_COMPILER_VERSION=\"$<CXX_COMPILER_VERSION>\"
)

target_link_libraries(
  hictk_internal_tools_build_metadata
  PRIVATE
    hictk::project_options
    hictk::project_warnings
)

set(
  HICTK_DEPENDENCIES
  Boost
  bshoshany-thread-pool
  CLI11
  concurrentqueue
  FastFloat
  fmt
  HDF5
  HighFive
  LibArchive
  libdeflate
  nlohmann_json
  phmap
  opentelemetry-cpp
  readerwriterqueue
  span-lite
  spdlog
  tomlplusplus
  zstd
)

foreach(PKG ${HICTK_DEPENDENCIES})
  find_package(${PKG} REQUIRED QUIET)
  string(TOUPPER ${PKG} PKG_UPPER)
  string(REGEX REPLACE "[^a-zA-Z0-9_]" _ PKG_NORMALIZED ${PKG_UPPER})
  target_compile_definitions(
    hictk_internal_tools_build_metadata
    PRIVATE
      "HICTK_${PKG_NORMALIZED}_VERSION=\"${${PKG}_VERSION}\""
  )
endforeach()

target_link_system_libraries(
    hictk_internal_tools_build_metadata
    PRIVATE
      HDF5::HDF5
      LibArchive::LibArchive
      "libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
      "zstd::libzstd_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>"
)

if(TARGET opentelemetry-cpp::opentelemetry-cpp)
  target_link_libraries(
    hictk_internal_tools_build_metadata
    PRIVATE
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::opentelemetry-cpp>"
  )
else()
  target_link_system_libraries(
    hictk_internal_tools_build_metadata
    PRIVATE
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::sdk>"
  )
endif()
