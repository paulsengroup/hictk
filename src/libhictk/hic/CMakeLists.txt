# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED)
find_package(FMT CONFIG REQUIRED)
find_package(libdeflate CONFIG REQUIRED)
find_package(phmap CONFIG REQUIRED)

if(HICTK_WITH_EIGEN)
  find_package(Eigen3 QUIET)
  if(NOT Eigen3_FOUND)
    find_package(
      Eigen3
      CONFIG
      QUIET
      REQUIRED)
  endif()
endif()

add_library(hic INTERFACE)
add_library(hictk::hic ALIAS hic)

target_sources(
  hic
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(hic INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                                         $<INSTALL_INTERFACE:include>)
target_link_libraries(
  hic
  INTERFACE hictk::balancing
            hictk::common
            hictk::chromosome
            hictk::reference)

target_link_system_libraries(
  hic
  INTERFACE
  $<$<BOOL:${HICTK_WITH_EIGEN}>:Eigen3::Eigen>
  fmt::fmt
  libdeflate::libdeflate_$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,shared,static>
  phmap
  std::filesystem)

include(GNUInstallDirs)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
