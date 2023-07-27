# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_subdirectory(balancing)
add_subdirectory(bin_table)
add_subdirectory(chromosome)
add_subdirectory(common)
add_subdirectory(cooler)
add_subdirectory(fmt)
add_subdirectory(genomic_interval)
add_subdirectory(hic)
add_subdirectory(numeric)
add_subdirectory(pixel)
add_subdirectory(reference)
add_subdirectory(transformers)
add_subdirectory(variant)

add_library(libhictk INTERFACE)
add_library(hictk::libhictk ALIAS libhictk)

target_link_libraries(
  libhictk
  INTERFACE hictk::balancing
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::cooler
            hictk::format
            hictk::genomic_interval
            hictk::hic
            hictk::numeric
            hictk::pixel
            hictk::reference
            hictk::transformers
            hictk::variant)

target_include_directories(
  libhictk
  INTERFACE
  INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>" "$<INSTALL_INTERFACE:include>")

install(
  TARGETS libhictk
          balancing
          bin_table
          chromosome
          common
          cooler
          formatting
          genomic_interval
          hic
          numeric
          pixel
          reference
          transformers
          variant
  EXPORT libhictk-targets
  FILE_SET HEADERS
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  PRIVATE_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  INCLUDES
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(
  EXPORT libhictk-targets
  FILE hictkTargets.cmake
  NAMESPACE hictk::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk/")

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/hictkConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/hictkConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk/")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/hictkConfig.cmake" DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk")
