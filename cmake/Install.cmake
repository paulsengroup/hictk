# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(
  targets
  balancing
  bin_table
  binary_buffer
  chromosome
  common/common
  common/default_delete_libdeflate
  common/default_delete_zstd
  common/genomic_units
  common/hash
  common/static_binary_buffer
  common/string
  common/tmpdir
  common/type_pretty_printer
  common/type_traits
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
  TARGETS
    hictk_internal_balancing
    hictk_internal_bin_table
    hictk_internal_binary_buffer
    hictk_internal_chromosome
    hictk_internal_common
    hictk_internal_common_default_delete_libdeflate
    hictk_internal_common_default_delete_zstd
    hictk_internal_common_genomic_units
    hictk_internal_common_hash
    hictk_internal_common_static_binary_buffer
    hictk_internal_common_string
    hictk_internal_common_tmpdir
    hictk_internal_common_type_pretty_printer
    hictk_internal_common_type_traits
    hictk_internal_common_version
    hictk_internal_cooler
    hictk_internal_expected_values_aggregator
    hictk_internal_file
    hictk_internal_filestream
    hictk_internal_formatting
    hictk_internal_genomic_interval
    hictk_internal_hic
    hictk_internal_numeric
    hictk_internal_pixel
    hictk_internal_reference
    hictk_internal_transformers
    hictk_internal_variant
  EXPORT libhictk-internal-targets
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

install(
  EXPORT libhictk-internal-targets
  COMPONENT Libraries
  FILE hictkInternalTargets.cmake
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
