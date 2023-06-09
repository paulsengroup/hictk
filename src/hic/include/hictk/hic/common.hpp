// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <tsl/ordered_map.h>

#include <cassert>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace hictk {

struct SerializedPixel {
  std::int64_t bin1_id{};
  std::int64_t bin2_id{};
  float count{};

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
  std::int64_t position{-1};
  std::int64_t size{-1};

  constexpr explicit operator bool() const noexcept { return size >= 0 && position >= 0; }
  constexpr bool operator<(const indexEntry &other) const noexcept {
    return position < other.position;
  }
  constexpr bool operator==(const indexEntry &other) const noexcept {
    return position == other.position && size == other.size;
  }
  constexpr bool operator!=(const indexEntry &other) const noexcept { return !(*this == other); }
};

}  // namespace hictk

namespace hictk {

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

}  // namespace hictk

template <>
struct fmt::formatter<hictk::NormalizationMethod> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::NormalizationMethod n, FormatContext &ctx)
      -> decltype(ctx.out()) {
    switch (n) {
      case hictk::NormalizationMethod::NONE:
        return fmt::format_to(ctx.out(), FMT_STRING("NONE"));
      case hictk::NormalizationMethod::VC:
        return fmt::format_to(ctx.out(), FMT_STRING("VC"));
      case hictk::NormalizationMethod::VC_SQRT:
        return fmt::format_to(ctx.out(), FMT_STRING("VC_SQRT"));
      case hictk::NormalizationMethod::KR:
        return fmt::format_to(ctx.out(), FMT_STRING("KR"));
      case hictk::NormalizationMethod::SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("SCALE"));
      case hictk::NormalizationMethod::INTER_VC:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_VC"));
      case hictk::NormalizationMethod::INTER_KR:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_KR"));
      case hictk::NormalizationMethod::INTER_SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("INTER_SCALE"));
      case hictk::NormalizationMethod::GW_VC:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_VC"));
      case hictk::NormalizationMethod::GW_KR:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_KR"));
      case hictk::NormalizationMethod::GW_SCALE:
        return fmt::format_to(ctx.out(), FMT_STRING("GW_SCALE"));
    }
    assert(false);
    std::abort();
  }
};

template <>
struct fmt::formatter<hictk::MatrixType> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::MatrixType t, FormatContext &ctx) -> decltype(ctx.out()) {
    switch (t) {
      case hictk::MatrixType::observed:
        return fmt::format_to(ctx.out(), FMT_STRING("observed"));
      case hictk::MatrixType::oe:
        return fmt::format_to(ctx.out(), FMT_STRING("oe"));
      case hictk::MatrixType::expected:
        return fmt::format_to(ctx.out(), FMT_STRING("expected"));
    }
    assert(false);
    std::abort();
  }
};

template <>
struct fmt::formatter<hictk::MatrixUnit> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::MatrixUnit u, FormatContext &ctx) -> decltype(ctx.out()) {
    switch (u) {
      case hictk::MatrixUnit::BP:
        return fmt::format_to(ctx.out(), FMT_STRING("BP"));
      case hictk::MatrixUnit::FRAG:
        return fmt::format_to(ctx.out(), FMT_STRING("FRAG"));
    }
    assert(false);
    std::abort();
  }
};

template <typename T>
using UniquePtrWithDeleter = std::unique_ptr<T, std::function<void(T *)>>;

struct GenomicCoordinates {
  std::string chrom;
  std::int32_t start;
  std::int32_t end;

  [[nodiscard]] inline static GenomicCoordinates fromString(std::string coord,
                                                            bool noChromName = false) {
    GenomicCoordinates gc{};

    const auto original_coord = coord;

    if (!noChromName) {
      auto pos = coord.find(':');
      if (pos == std::string::npos) {
        gc.chrom = coord;
        return gc;
      }

      gc.chrom = coord.substr(0, pos);
      coord = coord.substr(pos + 1);
    }

    auto pos = coord.find('-');
    if (pos == std::string::npos) {
      pos = coord.find(':');
    }
    if (pos == std::string::npos) {
      throw std::runtime_error(fmt::format(FMT_STRING("unable to parse coordinate \"{}\""), coord));
    }

    try {
      std::size_t tail{0};
      gc.start = std::stoi(coord.substr(0, pos));
      gc.end = std::stoi(coord.substr(pos + 1), &tail);
      if (gc.start >= gc.end) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid coordinate {}: start position >= end position"), coord));
      }
      if (gc.start < 0) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid coordinate {}: start position is negative"), coord));
      }
      coord = coord.substr(pos + 1);
      if (tail != coord.size()) {
        throw std::runtime_error(fmt::format(FMT_STRING("unable to parse \"{}\""), coord));
      }
    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(FMT_STRING("unable to parse coordinate \"{}\": {}"),
                                           original_coord, e.what()));
    }
    return gc;
  }
};
