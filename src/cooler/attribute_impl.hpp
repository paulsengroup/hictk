// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <highfive/H5Utility.hpp>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

#include "hictk/common.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_pretty_printer.hpp"

namespace hictk {

template <typename ParentObj>
inline bool Attribute::exists(ParentObj& h5obj, std::string_view key) {
  return h5obj.hasAttribute(std::string{key});
}

template <typename T, typename ParentObj>
inline void Attribute::write(ParentObj& h5obj, std::string_view key, const T& value,
                             bool overwrite_if_exists) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  const std::string key_{key};
  if (overwrite_if_exists && Attribute::exists(h5obj, key)) {
    h5obj.deleteAttribute(key_);
  }
  h5obj.createAttribute(key_, value);
}

template <typename T, typename ParentObj>
inline T Attribute::read(const ParentObj& h5obj, std::string_view key) {
  try {
    auto attrv = Attribute::read(h5obj, key, false);
    std::optional<T> buff{};
    std::visit(
        [&](auto& attr) {
          using Tin = remove_cvref_t<decltype(attr)>;
          using Tout = T;
          if constexpr (std::is_same_v<Tin, Tout>) {
            buff = attr;
          } else if constexpr (!std::is_same_v<Tin, std::monostate>) {
            buff = numeric_converter<Tin, Tout>(attr);
          }
        },
        attrv);

    assert(buff);

    return *buff;
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Unable to read attribute \"{}/{}\": {}"),
                                         h5obj.getPath(), key, e.what()));
  }
}

template <typename ParentObj>
inline auto Attribute::read(const ParentObj& h5obj, std::string_view key, bool missing_ok)
    -> AttributeVar {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT

  if (missing_ok && !Attribute::exists(h5obj, key)) {
    return std::monostate();
  }
  auto attr = read_variant(h5obj.getAttribute(std::string{key}));
  if (std::holds_alternative<std::monostate>(attr)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to read attribute \"{}\" from path \"{}\". Reason: "
                               "attribute exists but type is not supported"),
                    key, h5obj.getPath()));
  }
  return attr;
}

template <typename T, typename ParentObj>
inline std::vector<T> Attribute::read_vector(const ParentObj& h5obj, std::string_view key) {
  std::vector<T> buff;
  Attribute::read_vector(h5obj, key, buff);
  return buff;
}

template <typename T, typename ParentObj>
inline void Attribute::read_vector(const ParentObj& h5obj, std::string_view key,
                                   std::vector<T>& buff) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  try {
    h5obj.getAttribute(std::string{key}).read(buff);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Unable to read attribute \"{}/{}\": {}"),
                                         h5obj.getPath(), key, e.what()));
  }
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
template <std::size_t i>
inline auto Attribute::read_variant(const HighFive::Attribute& attr) -> AttributeVar {
  if constexpr (i < std::variant_size_v<AttributeVar>) {
    using T = std::variant_alternative_t<i, AttributeVar>;
    if (attr.getDataType() != HighFive::create_datatype<T>()) {
      return read_variant<i + 1>(attr);
    }
    T buff{};
    attr.read(buff);
    return buff;
  }
  return std::monostate();
}
DISABLE_WARNING_POP

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
template <typename T1, typename Tout, typename Tin>
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
inline Tout Attribute::numeric_converter(T1& buff) {
  static_assert(!std::is_same_v<Tin, std::monostate>);

  if constexpr (std::is_same_v<Tin, Tout>) {
    return buff;
  }

  if constexpr (std::is_floating_point_v<Tin> && std::is_floating_point_v<Tout>) {
    // Here we assume that errors due to casting e.g. float to double are acceptable
    return static_cast<Tout>(buff);
  }

  if constexpr (std::is_same_v<Tin, std::string> && std::is_arithmetic_v<Tout>) {
    // Try to convert a string attribute to the appropriate numeric type
    try {
      return internal::parse_numeric_or_throw<Tout>(buff);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected type {}, found std::string. An attempt to convert "
                                 "std::string to {} was made, but failed. Reason {}"),
                      internal::type_name<Tout>(), internal::type_name<Tin>(), e.what()));
    }
  }

  if constexpr (std::is_floating_point_v<Tin> && std::is_integral_v<Tout>) {
    // Only cast T2 to T1 only if conversion is lossless
    const auto lb = static_cast<Tin>(std::numeric_limits<Tout>::min());
    const auto ub = static_cast<Tin>(std::numeric_limits<Tout>::max());
    if (std::floor(buff) == buff && buff >= lb && buff <= ub) {
      return static_cast<Tout>(buff);
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Expected type {}, found {}. Unable to represent value {} as {} "
                               "without information loss"),
                    internal::type_name<Tout>(), internal::type_name<Tin>(), buff,
                    internal::type_name<Tout>()));
  }

  if constexpr (std::is_integral_v<Tin> && std::is_integral_v<Tout>) {
    // Cast integers without causing overflows
    if (buff < 0) {
      if (!std::is_unsigned_v<Tout>) {
        const auto lb = conditional_static_cast<std::int64_t>(std::numeric_limits<Tout>::min());
        const auto ub = conditional_static_cast<std::int64_t>(std::numeric_limits<Tout>::max());

        const auto buff_ = conditional_static_cast<std::int64_t>(buff);
        if (buff_ >= lb && buff_ <= ub) {
          return static_cast<Tout>(buff);
        }
      }
    } else {
      const auto ub = conditional_static_cast<std::uint64_t>(std::numeric_limits<Tout>::max());
      const auto buff_ = conditional_static_cast<std::uint64_t>(buff);
      if (buff_ <= ub) {
        return static_cast<Tout>(buff);
      }
    }
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Expected type {}, found {}. Unable to represent value {} as {} without overflowing"),
        internal::type_name<Tout>(), internal::type_name<Tin>(), buff,
        internal::type_name<Tout>()));
  }
  // No conversion was possible
  throw std::runtime_error(fmt::format(
      FMT_STRING(
          "Expected type {}, found {}. Unable to safely convert value {} of type {} to type {}"),
      internal::type_name<Tout>(), internal::type_name<Tin>(), buff, internal::type_name<Tin>(),
      internal::type_name<Tout>()));
}
DISABLE_WARNING_POP

}  // namespace hictk
