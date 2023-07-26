// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// This file was generated automatically by CMake.

#include <string_view>

namespace hictk::config::git {

// clang-format off
[[nodiscard]] constexpr bool state_available() noexcept { return true; }
[[nodiscard]] constexpr std::string_view head_sha1() noexcept { return "fc50422d004e96bfb16bb9b0ab6cb6f1246822a1"; }
[[nodiscard]] constexpr bool is_dirty() noexcept { return true; }
[[nodiscard]] constexpr std::string_view author_name() noexcept { return "Roberto Rossini"; }
[[nodiscard]] constexpr std::string_view author_email() noexcept { return "71787608+robomics@users.noreply.github.com"; }
[[nodiscard]] constexpr std::string_view commit_date() noexcept { return "2023-07-26 11:48:18 +0200"; }
[[nodiscard]] constexpr std::string_view commit_subject() noexcept { return "Avoid inclusing external header-only libs with add_subdirectory"; }
[[nodiscard]] constexpr std::string_view commit_body() noexcept { return "This makes it easier to exclude headers from install targets."; }
[[nodiscard]] constexpr std::string_view describe() noexcept { return "fc50422"; }
[[nodiscard]] constexpr std::string_view branch() noexcept { return "update-build-system"; }
[[nodiscard]] constexpr std::string_view tag() noexcept { return ""; }
// clang-format on

}  // namespace hictk::config::git
