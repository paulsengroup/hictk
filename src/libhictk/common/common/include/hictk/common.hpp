// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hictk {

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

[[nodiscard]] constexpr bool ndebug_defined() noexcept {
#ifdef NDEBUG
  return true;
#else
  return false;
#endif
}

[[nodiscard]] constexpr bool ndebug_not_defined() noexcept { return !ndebug_defined(); }

#if defined(__GNUC__)
#define HICTK_UNREACHABLE_CODE __builtin_unreachable()
#elif defined(_MSC_VER)
#define HICTK_UNREACHABLE_CODE __assume(0)
#else
#define HICTK_UNREACHABLE_CODE
#endif

[[noreturn]] inline void unreachable_code() {
  if constexpr (ndebug_not_defined()) {
    throw std::logic_error("Unreachable code");
  }
  HICTK_UNREACHABLE_CODE;
}

#if defined(__GNUC__)
#define HICTK_LIKELY(x)   __builtin_expect(!!(x), 1)  // NOLINT(*-simplify-boolean-expr)
#define HICTK_UNLIKELY(x) __builtin_expect(!!(x), 0)  // NOLINT(*-simplify-boolean-expr)
#else
#define HICTK_LIKELY(x)   x
#define HICTK_UNLIKELY(x) x
#endif

// NOLINTEND(cppcoreguidelines-macro-usage)

struct identity {
  template <typename T>
  [[nodiscard]] constexpr T &&operator()(T &&a) const noexcept {
    return std::forward<T>(a);
  }
  using is_transparent = void;
};

// to avoid useless casts (see https://github.com/nlohmann/json/issues/2893#issuecomment-889152324)
template <typename T, typename U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
  if constexpr (std::is_same_v<T, U>) {
    return value;
  } else {
    return static_cast<T>(value);
  }
}

}  // namespace hictk
