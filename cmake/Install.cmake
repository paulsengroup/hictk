# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(
  targets
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
)

include(GNUInstallDirs)
foreach(tgt ${targets})
  install(
    DIRECTORY
      "${PROJECT_SOURCE_DIR}/src/libhictk/${tgt}/include/hictk"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT Libraries
  )
endforeach()

install(
  TARGETS
    libhictk
    hictk_balancing
    hictk_bin_table
    hictk_binary_buffer
    hictk_chromosome
    hictk_common
    hictk_cooler
    hictk_expected_values_aggregator
    hictk_file
    hictk_filestream
    hictk_formatting
    hictk_genomic_interval
    hictk_hic
    hictk_numeric
    hictk_pixel
    hictk_reference
    hictk_transformers
    hictk_variant
  EXPORT libhictk-targets
  COMPONENT Libraries
  FILE_SET
  HEADERS
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  RUNTIME
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE
    DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  PUBLIC_HEADER
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  PRIVATE_HEADER
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(
  EXPORT libhictk-targets
  COMPONENT Libraries
  FILE hictkTargets.cmake
  NAMESPACE hictk::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk/"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/hictkConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/hictkConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk/"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/hictkConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/hictk"
  COMPONENT Libraries
)

install(
  FILES
    "${PROJECT_SOURCE_DIR}/LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/hictk/"
  COMPONENT Libraries
)
