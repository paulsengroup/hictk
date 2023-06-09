# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(FetchContent)

# cmake-format: off
FetchContent_Declare(
        _hictk_project_options
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/project_options-v0.29.0.tar.xz
        URL_HASH SHA512=5cef5c4aab91328d41d0c846db326075eebab0bdcd53919d0bde288eb3ff4544ec51116234af25cb93475bb1ed806c51cef757868a56ceafdcf08a553afb6bf5
)

# We're fetching phmap directly so that we can better control build options to avoid ASAN false positives
FetchContent_Declare(
        _hictk_phmap
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-hashmap-v1.3.11.tar.xz
        URL_HASH SHA512=f4e1a9388c2d996781482eebc2f02eed02a9b3cfe84707da97797f0aae5e83d6758ede5785f7f9fce502cf1ceef09590e084ee7d3362ca39bf5b143981af7486
)
# cmake-format: on

FetchContent_MakeAvailable(_hictk_project_options _hictk_phmap)
include_directories(SYSTEM ${_hictk_phmap_SOURCE_DIR})
include(${_hictk_project_options_SOURCE_DIR}/Index.cmake)
