// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// GCC to MSVC codes: https://github.com/srz-zumix/awesome-cpp-warning

// clang-format off

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

// Defines for MSVC
#ifdef _MSC_VER
    #define HICTK_DISABLE_WARNING_PUSH                      __pragma(warning(push))
    #define HICTK_DISABLE_WARNING_POP                       __pragma(warning(pop))
    #define DISABLE_WARNING(warningNumber)            __pragma(warning(disable : warningNumber))

    #define HICTK_DISABLE_WARNING_BOOL_COMPARE              DISABLE_WARNING(4806)
    #define HICTK_DISABLE_WARNING_CONVERSION                DISABLE_WARNING(C4244)
    #define HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS   DISABLE_WARNING(4996)
    #define HICTK_DISABLE_WARNING_MAYBE_UNINITIALIZED
    #define HICTK_DISABLE_WARNING_NULL_DEREF
    #define HICTK_DISABLE_WARNING_USELESS_CAST
    #define HICTK_DISABLE_WARNING_SIGN_COMPARE
    #define HICTK_DISABLE_WARNING_UNREACHABLE_CODE          DISABLE_WARNING(4702)
#endif

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define DO_PRAGMA(X)                              _Pragma(#X)
    #define HICTK_DISABLE_WARNING_PUSH                      DO_PRAGMA(GCC diagnostic push)
    #define HICTK_DISABLE_WARNING_POP                       DO_PRAGMA(GCC diagnostic pop)
    #define DISABLE_WARNING(warningName)              DO_PRAGMA(GCC diagnostic ignored warningName)

    #define HICTK_DISABLE_WARNING_CONVERSION                DISABLE_WARNING("-Wconversion")
    #define HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS   DISABLE_WARNING("-Wdeprecated-declarations")
    #define HICTK_DISABLE_WARNING_NULL_DEREF                DISABLE_WARNING("-Wnull-dereference")
    #define HICTK_DISABLE_WARNING_SIGN_COMPARE              DISABLE_WARNING("-Wsign-compare")
    #define HICTK_DISABLE_WARNING_UNREACHABLE_CODE
#endif

// Defines for GCC only
#if defined(__GNUC__) && !defined(__clang__)
    #define HICTK_DISABLE_WARNING_BOOL_COMPARE              DISABLE_WARNING("-Wbool-compare")
    #define HICTK_DISABLE_WARNING_MAYBE_UNINITIALIZED       DISABLE_WARNING("-Wmaybe-uninitialized")
    #define HICTK_DISABLE_WARNING_USELESS_CAST              DISABLE_WARNING("-Wuseless-cast")
#endif

// Defines for Clang only
#ifdef __clang__
    #define HICTK_DISABLE_WARNING_BOOL_COMPARE              DISABLE_WARNING("-Wtautological-constant-out-of-range-compare")
    #define HICTK_DISABLE_WARNING_MAYBE_UNINITIALIZED
    #define HICTK_DISABLE_WARNING_USELESS_CAST
#endif

// Defines for unknown/unsupported compilers
#if !defined(_MSC_VER) && !defined(__GNUC__) && !defined(__clang__)
  #define DISABLE_WARNING
  #define HICTK_DISABLE_WARNING_PUSH
  #define HICTK_DISABLE_WARNING_POP

  #define HICTK_DISABLE_WARNING_BOOL_COMPARE
  #define HICTK_DISABLE_WARNING_CONVERSION
  #define HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
  #define HICTK_DISABLE_WARNING_MAYBE_UNINITIALIZED
  #define HICTK_DISABLE_WARNING_NULL_DEREF
  #define HICTK_DISABLE_WARNING_USELESS_CAST
  #define HICTK_DISABLE_WARNING_SIGN_COMPARE
  #define HICTK_DISABLE_WARNING_UNREACHABLE_CODE
#endif

// NOLINTEND(cppcoreguidelines-macro-usage)

// clang-format on
