// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <string_view>
#include <type_traits>

namespace hictk {

template <typename T>
struct remove_cvref {
  using type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

template <typename T>
struct is_string
    : public std::disjunction<std::is_same<char *, typename std::decay_t<T>>,
                              std::is_same<const char *, typename std::decay_t<T>>,
                              std::is_same<std::string, typename std::decay_t<T>>,
                              std::is_same<std::string_view, typename std::decay_t<T>>> {};

template <typename T>
constexpr bool is_string_v = is_string<T>::value;

template <typename T, typename Enabler = void>
struct is_map : std::false_type {};

template <typename T>
struct is_map<T, std::void_t<typename T::mapped_type>> : std::true_type {};

template <typename T>
constexpr bool is_map_v = is_map<T>::value;

template <typename Operation, typename Operand>
struct is_unary_operation : public std::is_invocable<Operation, Operand> {};

template <typename Operation, typename Operand>
constexpr bool is_unary_operation_v = is_unary_operation<Operation, Operand>::value;

namespace internal {
// Adapted from https://stackoverflow.com/a/29634934
// clang-format off
    template <typename T>
    auto is_iterable_impl(int)
    -> decltype(std::begin(std::declval<T &>()) != std::end(std::declval<T &>()),  // begin/end and operator !=
    void(),                                                            // Handle evil operator ,
    ++std::declval<decltype(std::begin(std::declval<T &>())) &>(),     // operator ++
    void(*std::begin(std::declval<T &>())),                            // operator*
    std::true_type{});
// clang-format on

template <typename T>
std::false_type is_iterable_impl(...);
}  // namespace internal

template <typename T, typename = void>
struct is_iterable : std::false_type {};

// clang-format off
template <typename T>
struct is_iterable<T, std::void_t<decltype(internal::is_iterable_impl<T>(0))>>
: std::true_type {};
// clang-format on

template <typename T>
constexpr bool is_iterable_v = is_iterable<T>::value;

template <typename Operation, typename It, typename = void>
constexpr bool is_unary_operation_on_iterator = false;

template <typename Operation, typename It>
constexpr bool is_unary_operation_on_iterator<
    Operation, It,
    std::void_t<std::disjunction<is_iterable<It>, std::is_pointer<It>>,
                is_unary_operation<Operation, decltype(*std::declval<It>())>>> = true;

}  // namespace hictk
