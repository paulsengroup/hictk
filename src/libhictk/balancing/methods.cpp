// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/methods.hpp"

#include <fmt/format.h>

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace hictk::balancing {

Method::Method(std::string_view name) : _name(std::string{name}) {
  if (_name.empty()) {
    throw std::runtime_error("weight dataset name is empty");
  }
}
bool operator==(const Method& a, const Method& b) noexcept { return a._name == b._name; }
bool operator==(const Method& a, std::string_view b) noexcept { return a._name == b; }
bool operator==(std::string_view a, const Method& b) noexcept { return b == a; }
bool operator!=(const Method& a, const Method& b) noexcept { return !(a == b); }
bool operator!=(const Method& a, std::string_view b) noexcept { return !(a == b); }
bool operator!=(std::string_view a, const Method& b) noexcept { return !(a == b); }

bool operator<(const Method& a, const Method& b) noexcept { return a._name < b._name; }
bool operator<(const Method& a, std::string_view b) noexcept { return a._name < b; }
bool operator<(std::string_view a, const Method& b) noexcept { return a < b._name; }

std::string_view Method::to_string() const noexcept { return _name; }

// NOLINTBEGIN(bugprone-exception-escape)
const Method& Method::NONE() noexcept {
  static const Method m{"NONE"};
  return m;
}
const Method& Method::VC() noexcept {
  static const Method m{"VC"};
  return m;
}
const Method& Method::VC_SQRT() noexcept {
  static const Method m{"VC_SQRT"};
  return m;
}
const Method& Method::KR() noexcept {
  static const Method m{"KR"};
  return m;
}
const Method& Method::SCALE() noexcept {
  static const Method m{"SCALE"};
  return m;
}
const Method& Method::ICE() noexcept {
  static const Method m{"ICE"};
  return m;
}
const Method& Method::INTER_VC() noexcept {
  static const Method m{"INTER_VC"};
  return m;
}
const Method& Method::INTER_KR() noexcept {
  static const Method m{"INTER_KR"};
  return m;
}
const Method& Method::INTER_SCALE() noexcept {
  static const Method m{"INTER_SCALE"};
  return m;
}
const Method& Method::INTER_ICE() noexcept {
  static const Method m{"INTER_ICE"};
  return m;
}
const Method& Method::GW_VC() noexcept {
  static const Method m{"GW_VC"};
  return m;
}
const Method& Method::GW_KR() noexcept {
  static const Method m{"GW_KR"};
  return m;
}
const Method& Method::GW_SCALE() noexcept {
  static const Method m{"GW_SCALE"};
  return m;
}
const Method& Method::GW_ICE() noexcept {
  static const Method m{"GW_ICE"};
  return m;
}
// NOLINTEND(bugprone-exception-escape)
}  // namespace hictk::balancing

std::size_t std::hash<hictk::balancing::Method>::operator()(
    const hictk::balancing::Method& m) const noexcept {
  return std::hash<std::string_view>{}(m.to_string());
}
