// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string_view>

namespace hictk::internal {
// Source: https://stackoverflow.com/a/56600402

template <typename T>
[[nodiscard]] constexpr std::string_view type_name();

template <>
[[nodiscard]] constexpr std::string_view type_name<void>() {
  return "void";
}

namespace detail {

using type_name_prober = void;

template <typename T>
constexpr std::string_view wrapped_type_name() {
#ifdef __clang__
  return __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
  return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
  return __FUNCSIG__;
#else
#error "Unsupported compiler"
#endif
}

constexpr std::size_t wrapped_type_name_prefix_length() {
  return wrapped_type_name<type_name_prober>().find(type_name<type_name_prober>());
}

constexpr std::size_t wrapped_type_name_suffix_length() {
  return wrapped_type_name<type_name_prober>().length() - wrapped_type_name_prefix_length() -
         type_name<type_name_prober>().length();
}

}  // namespace detail

template <typename T>
constexpr std::string_view type_name() {
  constexpr auto wrapped_name = detail::wrapped_type_name<T>();
  constexpr auto prefix_length = detail::wrapped_type_name_prefix_length();
  constexpr auto suffix_length = detail::wrapped_type_name_suffix_length();
  constexpr auto type_name_length = wrapped_name.length() - prefix_length - suffix_length;
  return wrapped_name.substr(prefix_length, type_name_length);
}
}  // namespace hictk::internal
