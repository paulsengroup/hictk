// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include "hictk/common.hpp"

namespace hictk::hic {

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
  std::int64_t position{-1};  // NOLINT
  std::int64_t size{-1};      // NOLINT

  constexpr explicit operator bool() const noexcept { return size >= 0 && position >= 0; }
  constexpr bool operator<(const indexEntry &other) const noexcept {
    return position < other.position;
  }
  constexpr bool operator==(const indexEntry &other) const noexcept {
    return position == other.position && size == other.size;
  }
  constexpr bool operator!=(const indexEntry &other) const noexcept { return !(*this == other); }
};

}  // namespace hictk::hic

namespace hictk::hic {

enum class MatrixType { observed, oe, expected };
enum class MatrixUnit { BP, FRAG };

[[nodiscard]] inline MatrixType ParseMatrixTypeStr(const std::string &s) {
  if (s == "observed") {
    return MatrixType::observed;
  }
  if (s == "oe") {
    return MatrixType::oe;
  }
  if (s == "expected") {
    return MatrixType::expected;
  }

  throw std::runtime_error("Invalid matrix type \"" + s + "\"");
}

[[nodiscard]] inline MatrixUnit ParseUnitStr(const std::string &s) {
  if (s == "BP") {
    return MatrixUnit::BP;
  }
  if (s == "FRAG") {
    return MatrixUnit::FRAG;
  }

  throw std::runtime_error("Invalid unit \"" + s + "\"");
}

}  // namespace hictk::hic

template <>
struct fmt::formatter<hictk::hic::MatrixType> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::hic::MatrixType t, FormatContext &ctx) -> decltype(ctx.out()) {
    switch (t) {
      case hictk::hic::MatrixType::observed:
        return fmt::format_to(ctx.out(), FMT_STRING("observed"));
      case hictk::hic::MatrixType::oe:
        return fmt::format_to(ctx.out(), FMT_STRING("oe"));
      case hictk::hic::MatrixType::expected:
        return fmt::format_to(ctx.out(), FMT_STRING("expected"));
    }
    HICTK_UNREACHABLE_CODE;
  }
};

template <>
struct fmt::formatter<hictk::hic::MatrixUnit> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::hic::MatrixUnit u, FormatContext &ctx) -> decltype(ctx.out()) {
    switch (u) {
      case hictk::hic::MatrixUnit::BP:
        return fmt::format_to(ctx.out(), FMT_STRING("BP"));
      case hictk::hic::MatrixUnit::FRAG:
        return fmt::format_to(ctx.out(), FMT_STRING("FRAG"));
    }
    HICTK_UNREACHABLE_CODE;
  }
};

template <typename T>
using UniquePtrWithDeleter = std::unique_ptr<T, std::function<void(T *)>>;
