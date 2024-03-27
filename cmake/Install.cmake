# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(targets
    balancing;bin_table;binary_buffer;chromosome;common;cooler;expected_values_aggregator;file;filestream;formatting;genomic_interval;hic;numeric;pixel;reference;transformers;variant
)

include(GNUInstallDirs)
foreach(tgt ${targets})
  install(DIRECTORY "${PROJECT_SOURCE_DIR}/src/libhictk/${tgt}/include/hictk" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endforeach()

install(
  TARGETS libhictk
          balancing
          bin_table
          binary_buffer
          chromosome
          common
          cooler
          expected_values_aggregator
          file
          filestream
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
