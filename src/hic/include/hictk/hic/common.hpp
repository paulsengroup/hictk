// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "hictk/common.hpp"

namespace hictk::hic {

struct SerializedPixel {
  std::int64_t bin1_id{};  // NOLINT
  std::int64_t bin2_id{};  // NOLINT
  float count{};           // NOLINT

  constexpr bool operator<(const SerializedPixel &other) const noexcept {
    if (bin1_id == other.bin1_id) {
      return bin2_id < other.bin2_id;
    }
    return bin1_id < other.bin1_id;
  }
  constexpr bool operator==(const SerializedPixel &other) const noexcept {
    return bin1_id == other.bin1_id && bin2_id == other.bin2_id && count == other.count;
  }
  constexpr bool operator!=(const SerializedPixel &other) const noexcept {
    return !(*this == other);
  }
};

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

enum class NormalizationMethod {
  NONE,
  VC,
  VC_SQRT,
  KR,
  SCALE,
  INTER_VC,
  INTER_KR,
  INTER_SCALE,
  GW_VC,
  GW_KR,
  GW_SCALE
};
enum class MatrixType { observed, oe, expected };
enum class MatrixUnit { BP, FRAG };

[[nodiscard]] inline NormalizationMethod ParseNormStr(const std::string &s) {
  if (s == "NONE") {
    return NormalizationMethod::NONE;
  }
  if (s == "VC") {
    return NormalizationMethod::VC;
  }

  if (s == "VC_SQRT") {
    return NormalizationMethod::VC_SQRT;
  }

  if (s == "KR") {
    return NormalizationMethod::KR;
  }

  if (s == "SCALE") {
    return NormalizationMethod::SCALE;
  }

  if (s == "INTER_VC") {
    return NormalizationMethod::INTER_VC;
  }

  if (s == "INTER_KR") {
    return NormalizationMethod::INTER_KR;
  }

  if (s == "INTER_SCALE") {
    return NormalizationMethod::INTER_SCALE;
  }

  if (s == "GW_VC") {
    return NormalizationMethod::GW_VC;
  }

  if (s == "GW_KR") {
    return NormalizationMethod::GW_KR;
  }

  if (s == "GW_SCALE") {
    return NormalizationMethod::GW_SCALE;
  }

  throw std::runtime_error("Invalid normalization \"" + s + "\"");
}

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
struct fmt::formatter<hictk::hic::NormalizationMethod> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::hic::NormalizationMethod n, FormatContext &ctx)
      -> decltype(ctx.out()) {
    using NM = hictk::hic::NormalizationMethod;
    switch (n) {
      case NM::NONE:
        return fmt::format_to(ctx.out(), FMT_STRING("NONE"));
      case NM::VC:
        return fmt::format_to(ctx.out(), FMT_STRING("VC"));
      case NM::VC_SQRT:
        return fmt::format_to(ctx.out(), FMT_STRING("VC_SQRT"));
      case NM::KR:
        return fmt::format_to(ctx.out(), FMT_STRING("KR"));
      case NM::SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("SCALE"));
      case NM::INTER_VC:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_VC"));
      case NM::INTER_KR:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_KR"));
      case NM::INTER_SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_SCALE"));
      case NM::GW_VC:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_VC"));
      case NM::GW_KR:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_KR"));
      case NM::GW_SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_SCALE"));
    }
    HICTK_UNREACHABLE_CODE;
  }
};

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
