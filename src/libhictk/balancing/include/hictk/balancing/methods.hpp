// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <string>
#include <string_view>

namespace hictk::balancing {

class Method {
  std::string _name{};

 public:
  Method() = default;
  explicit Method(std::string_view name);
  friend bool operator==(const Method &a, const Method &b) noexcept;
  friend bool operator==(const Method &a, std::string_view b) noexcept;
  friend bool operator==(std::string_view a, const Method &b) noexcept;
  friend bool operator!=(const Method &a, const Method &b) noexcept;
  friend bool operator!=(const Method &a, std::string_view b) noexcept;
  friend bool operator!=(std::string_view a, const Method &b) noexcept;

  friend bool operator<(const Method &a, const Method &b) noexcept;
  friend bool operator<(const Method &a, std::string_view b) noexcept;
  friend bool operator<(std::string_view a, const Method &b) noexcept;

  [[nodiscard]] std::string_view to_string() const noexcept;

  static const Method &NONE() noexcept;
  static const Method &VC() noexcept;
  static const Method &VC_SQRT() noexcept;
  static const Method &KR() noexcept;
  static const Method &SCALE() noexcept;
  static const Method &ICE() noexcept;
  static const Method &INTER_VC() noexcept;
  static const Method &INTER_KR() noexcept;
  static const Method &INTER_SCALE() noexcept;
  static const Method &INTER_ICE() noexcept;
  static const Method &GW_VC() noexcept;
  static const Method &GW_KR() noexcept;
  static const Method &GW_SCALE() noexcept;
  static const Method &GW_ICE() noexcept;
};
}  // namespace hictk::balancing

template <>
struct fmt::formatter<hictk::balancing::Method> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::balancing::Method &n, FormatContext &ctx) -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("{}"), n.to_string());
  }
};

template <>
struct std::hash<hictk::balancing::Method> {
  std::size_t operator()(const hictk::balancing::Method &m) const noexcept;
};
