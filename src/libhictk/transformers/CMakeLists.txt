# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(HICTK_WITH_ARROW)
  find_package(Arrow REQUIRED QUIET)
endif()

add_library(transformers INTERFACE)
add_library(hictk::transformers ALIAS transformers)

target_sources(
  transformers
  INTERFACE FILE_SET
            HEADERS
            BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(transformers INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                  "$<INSTALL_INTERFACE:include>")
target_link_libraries(transformers INTERFACE hictk::cooler hictk::hic)


if (BUILD_SHARED_LIBS)
  target_link_system_libraries(transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_shared>")
else()
  target_link_system_libraries(transformers INTERFACE "$<$<BOOL:${HICTK_WITH_ARROW}>:Arrow::arrow_static>")
endif()
