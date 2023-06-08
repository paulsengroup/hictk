// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// This file was generated automatically by CMake.

#include <string_view>

namespace hictk::config::git {

// clang-format off
[[nodiscard]] constexpr bool state_available() noexcept { return true; }
[[nodiscard]] constexpr std::string_view head_sha1() noexcept { return "18bc8f5f99d5763f33bcc054dc9ee0dc340d902d"; }
[[nodiscard]] constexpr bool is_dirty() noexcept { return true; }
[[nodiscard]] constexpr std::string_view author_name() noexcept { return "Roberto Rossini"; }
[[nodiscard]] constexpr std::string_view author_email() noexcept { return "71787608+robomics@users.noreply.github.com"; }
[[nodiscard]] constexpr std::string_view commit_date() noexcept { return "2023-06-08 12:54:13 +0200"; }
[[nodiscard]] constexpr std::string_view commit_subject() noexcept { return "Add code from hicxx"; }
[[nodiscard]] constexpr std::string_view commit_body() noexcept { return "https://github.com/robomics/hicxx"; }
[[nodiscard]] constexpr std::string_view describe() noexcept { return "18bc8f5"; }
[[nodiscard]] constexpr std::string_view branch() noexcept { return "main"; }
[[nodiscard]] constexpr std::string_view tag() noexcept { return ""; }
// clang-format on

}  // namespace hictk::config::git
