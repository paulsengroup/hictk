# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git" AND HICTK_ENABLE_GIT_VERSION_TRACKING)
  message(
    WARNING
    "-- Unable to find .git/ under \"${PROJECT_SOURCE_DIR}\". Setting -DHICTK_ENABLE_GIT_VERSION_TRACKING=OFF"
  )
  set(HICTK_ENABLE_GIT_VERSION_TRACKING OFF)
endif()

function(ConfigureVersioning input_config_folder output_config_folder)
  set(PRE_CONFIGURE_FILE "${input_config_folder}/git.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/git.hpp")

  if(HICTK_ENABLE_GIT_VERSION_TRACKING)
    include(FetchContent)
    FetchContent_Declare(
      _hictk_cmake-git-version-tracking
      URL
        "${PROJECT_SOURCE_DIR}/external/cmake-git-version-tracking.20250308.tar.xz"
      URL_HASH SHA256=524d80b76125ef1edf6156cd897171e587779f4521c6f3a80a0420336724ab1c
    )
    FetchContent_MakeAvailable(_hictk_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    include("${_hictk_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake")
  else()
    # Add dummy target
    add_custom_target(_hictk_check_git)

    if(NOT DEFINED HICTK_GIT_RETRIEVED_STATE)
      set(HICTK_GIT_RETRIEVED_STATE false)
    endif()
    if(NOT DEFINED HICTK_GIT_HEAD_SHA1)
      set(HICTK_GIT_HEAD_SHA1 "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_IS_DIRTY)
      set(HICTK_GIT_IS_DIRTY false)
    endif()
    if(NOT DEFINED HICTK_GIT_AUTHOR_NAME)
      set(HICTK_GIT_AUTHOR_NAME "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_AUTHOR_EMAIL)
      set(HICTK_GIT_AUTHOR_EMAIL "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_COMMIT_DATE_ISO8601)
      set(HICTK_GIT_COMMIT_DATE_ISO8601 "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_COMMIT_SUBJECT)
      set(HICTK_GIT_COMMIT_SUBJECT "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_COMMIT_BODY)
      set(HICTK_GIT_COMMIT_BODY "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_DESCRIBE)
      set(HICTK_GIT_DESCRIBE "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_BRANCH)
      set(HICTK_GIT_BRANCH "unknown")
    endif()
    if(NOT DEFINED HICTK_GIT_TAG)
      set(HICTK_GIT_TAG "unknown")
    endif()

    if(NOT WIN32)
      file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
    endif()
    configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
  endif()

  set(PRE_CONFIGURE_FILE "${input_config_folder}/version.hpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/version.hpp")

  if(NOT WIN32)
    file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
  endif()
  configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
endfunction()
