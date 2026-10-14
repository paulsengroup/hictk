// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

#include <fmt/format.h>

#include <cassert>
#include <exception>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <variant>

#include "hictk/cooler/attribute.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

constexpr RootGroup::operator bool() const noexcept { return _group.has_value(); }
constexpr bool RootGroup::operator!() const noexcept { return !bool(*this); }

constexpr HighFive::Group &RootGroup::operator()() {
  assert(_group.has_value());
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return _group.value();
}
constexpr const HighFive::Group &RootGroup::operator()() const {
  assert(_group.has_value());
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return _group.value();
}

template <typename N>
inline N RootGroup::read_attribute(std::string_view key) const {
  assert(_group.has_value());
  N buff;
  read_attribute(key, buff);
  return buff;
}

template <typename N>
inline void RootGroup::read_attribute(std::string_view key, N &buff) const {
  assert(_group.has_value());
  if (!try_read_attribute(key, buff)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Failed to read attribute \"{}\" from path \"{}\". Reason: attribute does not exist"),
        key, _group->getPath()));
  }
}

template <typename N>
inline bool RootGroup::try_read_attribute(std::string_view key, N &buff) const {
  assert(_group.has_value());
  // NOLINTBEGIN(bugprone-unchecked-optional-access)
  if (!Attribute::exists(*_group, key)) {
    return false;
  }

  try {
    using T = remove_cvref_t<decltype(*buff)>;
    buff = Attribute::read<T>(*_group, key);
    return true;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"), key,
                    _group->getPath(), e.what()));
  }
  // NOLINTEND(bugprone-unchecked-optional-access)
}

template <typename N>
inline std::optional<N> RootGroup::try_read_attribute(std::string_view key) const {
  N buff;
  if (try_read_attribute(key, buff)) {
    return buff;
  }

  return {};
}

template <typename N>
inline std::optional<N> RootGroup::try_read_sum_attribute(std::string_view key) const {
  N buff;
  if (try_read_attribute(key, buff)) {
    return buff;
  }

  return {};
}

template <typename N>
inline bool RootGroup::try_read_sum_attribute(std::string_view key, N &buff) const {
  assert(_group.has_value());
  // NOLINTBEGIN(bugprone-unchecked-optional-access)
  if (!Attribute::exists(*_group, key)) {
    return false;
  }

  try {
    auto sumv = Attribute::read(*_group, key);
    assert(!sumv.valueless_by_exception());
    const auto ok = std::visit(
        [&](auto sum) {  // NOLINT(performance-unnecessary-value-param)
          using T = remove_cvref_t<decltype(sum)>;
          if constexpr (std::is_integral_v<T>) {
            buff = conditional_static_cast<std::int64_t>(sum);
            return true;
          }
          if constexpr (std::is_floating_point_v<T>) {
            buff = conditional_static_cast<double>(sum);
            return true;
          }
          return false;
        },
        sumv);

    if (!ok) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("attribute \"{}{}\" does not have a numeric type"), _group->getPath(), key));
    }

    return ok;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"), key,
                    _group->getPath(), e.what()));
  }
  // NOLINTEND(bugprone-unchecked-optional-access)
}

constexpr HighFive::Group &Group::operator()() {
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return _group.value();
}
constexpr const HighFive::Group &Group::operator()() const {
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return _group.value();
}

}  // namespace hictk::cooler
