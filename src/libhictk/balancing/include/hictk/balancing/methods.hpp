// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <functional>
#include <string>
#include <string_view>
#include <utility>

namespace hictk::balancing {

class Method {
  std::string _name{};

 public:
  Method() = default;
  explicit Method(std::string_view name) : _name(std::string{name}) {
    if (_name.empty()) {
      throw std::runtime_error("weight dataset name is empty");
    }
  }

  [[nodiscard]] friend bool operator==(const Method &a, const Method &b) {
    return a._name == b._name;
  }
  [[nodiscard]] friend bool operator==(const Method &a, std::string_view b) { return a._name == b; }
  [[nodiscard]] friend bool operator==(std::string_view a, const Method &b) { return b == a; }
  [[nodiscard]] friend bool operator!=(const Method &a, const Method &b) { return !(a == b); }
  [[nodiscard]] friend bool operator!=(const Method &a, std::string_view b) { return !(a == b); }
  [[nodiscard]] friend bool operator!=(std::string_view a, const Method &b) { return !(a == b); }

  [[nodiscard]] std::string_view to_string() const noexcept { return _name; }

  static Method NONE() { return Method{"NONE"}; }
  static Method VC() { return Method{"VC"}; }
  static Method VC_SQRT() { return Method{"VC_SQRT"}; }
  static Method KR() { return Method{"KR"}; }
  static Method SCALE() { return Method{"SCALE"}; }
  static Method ICE() { return Method{"ICE"}; }
  static Method INTER_VC() { return Method{"INTER_VC"}; }
  static Method INTER_KR() { return Method{"INTER_KR"}; }
  static Method INTER_SCALE() { return Method{"INTER_SCALE"}; }
  static Method INTER_ICE() { return Method{"INTER_ICE"}; }
  static Method GW_VC() { return Method{"GW_VC"}; }
  static Method GW_KR() { return Method{"GW_KR"}; }
  static Method GW_SCALE() { return Method{"GW_SCALE"}; }
  static Method GW_ICE() { return Method{"GW_ICE"}; }
};
}  // namespace hictk::balancing

template <>
struct fmt::formatter<hictk::balancing::Method> {
  static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
      throw fmt::format_error("invalid format");
    }
    return ctx.end();
  }

  template <class FormatContext>
  static auto format(const hictk::balancing::Method n, FormatContext &ctx) -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("{}"), n.to_string());
  }
};

template <>
struct std::hash<hictk::balancing::Method> {
  inline std::size_t operator()(const hictk::balancing::Method &m) const noexcept {
    return std::hash<std::string_view>{}(m.to_string());
  }
};
