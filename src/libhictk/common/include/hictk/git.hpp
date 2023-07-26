// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// This file was generated automatically by CMake.

#include <string_view>

namespace hictk::config::git {

// clang-format off
[[nodiscard]] constexpr bool state_available() noexcept { return true; }
[[nodiscard]] constexpr std::string_view head_sha1() noexcept { return "53339ba40a2f2d7e054c269c6c45c2bd85567d31"; }
[[nodiscard]] constexpr bool is_dirty() noexcept { return true; }
[[nodiscard]] constexpr std::string_view author_name() noexcept { return "Roberto Rossini"; }
[[nodiscard]] constexpr std::string_view author_email() noexcept { return "71787608+robomics@users.noreply.github.com"; }
[[nodiscard]] constexpr std::string_view commit_date() noexcept { return "2023-07-26 13:12:49 +0200"; }
[[nodiscard]] constexpr std::string_view commit_subject() noexcept { return "Let Conan manage external deps"; }
[[nodiscard]] constexpr std::string_view commit_body() noexcept { return ""; }
[[nodiscard]] constexpr std::string_view describe() noexcept { return "53339ba"; }
[[nodiscard]] constexpr std::string_view branch() noexcept { return "update-build-system"; }
[[nodiscard]] constexpr std::string_view tag() noexcept { return ""; }
// clang-format on

}  // namespace hictk::config::git
