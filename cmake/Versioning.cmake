# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set(HICTK_PROJECT_VERSION_MAJOR 0)
set(HICTK_PROJECT_VERSION_MINOR 0)
set(HICTK_PROJECT_VERSION_PATCH 1)
set(HICTK_PROJECT_VERSION_SUFFIX "")

option(HICTK_ENABLE_GIT_VERSION_TRACKING "Retrieve project version and metadata from git" ON)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git" AND HICTK_ENABLE_GIT_VERSION_TRACKING)
  message(
    WARNING
      "-- Unable to find .git/ under \"${CMAKE_CURRENT_SOURCE_DIR}\". Setting -DHICTK_ENABLE_GIT_VERSION_TRACKING=OFF")
  set(HICTK_ENABLE_GIT_VERSION_TRACKING
      OFF
      CACHE FORCE)
endif()

function(ConfigureVersioning input_config_folder output_config_folder)
  set(PRE_CONFIGURE_FILE "${input_config_folder}/git.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/git.hpp")
  file(TOUCH ${POST_CONFIGURE_FILE})

  if(HICTK_ENABLE_GIT_VERSION_TRACKING)
    # cmake-format: off
    FetchContent_Declare(
            _hictk_cmake-git-version-tracking
            URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cmake-git-version-tracking.20221027.tar.xz
            URL_HASH SHA512=aa339b9fb5f147c5d88341ed4abf85e88b78173714c3a994860aa81f3b558b674829c17a1567d04642ae9df3ce6ed0e88001ba143bb152e91ed7e75bc607a86b
    )
    # cmake-format: on
    FetchContent_MakeAvailable(_hictk_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    include(${_hictk_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake)
  else()
    # Add dummy target
    add_custom_target(check_git)

    if(NOT DEFINED GIT_RETRIEVED_STATE)
      set(GIT_RETRIEVED_STATE false)
    endif()
    if(NOT DEFINED GIT_HEAD_SHA1)
      set(GIT_HEAD_SHA1 "unknown")
    endif()
    if(NOT DEFINED GIT_IS_DIRTY)
      set(GIT_IS_DIRTY false)
    endif()
    if(NOT DEFINED GIT_AUTHOR_NAME)
      set(GIT_AUTHOR_NAME "unknown")
    endif()
    if(NOT DEFINED GIT_AUTHOR_EMAIL)
      set(GIT_AUTHOR_EMAIL "unknown")
    endif()
    if(NOT DEFINED GIT_COMMIT_DATE_ISO8601)
      set(GIT_COMMIT_DATE_ISO8601 "unknown")
    endif()
    if(NOT DEFINED GIT_COMMIT_SUBJECT)
      set(GIT_COMMIT_SUBJECT "unknown")
    endif()
    if(NOT DEFINED GIT_COMMIT_BODY)
      set(GIT_COMMIT_BODY "unknown")
    endif()
    if(NOT DEFINED GIT_DESCRIBE)
      set(GIT_DESCRIBE "unknown")
    endif()
    if(NOT DEFINED GIT_BRANCH)
      set(GIT_BRANCH "unknown")
    endif()
    if(NOT DEFINED GIT_TAG)
      set(GIT_TAG "unknown")
    endif()

    configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE} @ONLY)
  endif()

  set(PRE_CONFIGURE_FILE "${input_config_folder}/version.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/version.hpp")

  file(TOUCH ${POST_CONFIGURE_FILE})
  configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE} @ONLY)
endfunction()

configureversioning("${CMAKE_CURRENT_SOURCE_DIR}/src/config" "${CMAKE_CURRENT_SOURCE_DIR}/src/common/include/hictk")
