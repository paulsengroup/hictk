// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/
// GCC to MSVC codes: https://github.com/srz-zumix/awesome-cpp-warning

// clang-format off

// Defines for MSVC
#ifdef _MSC_VER
    #define DISABLE_WARNING_PUSH                    __pragma(warning(push))
    #define DISABLE_WARNING_POP                     __pragma(warning(pop))
    #define DISABLE_WARNING(warningNumber)          __pragma(warning(disable : warningNumber))

    #define DISABLE_WARNING_DEPRECATED_DECLARATIONS
    #define DISABLE_WARNING_NULL_DEREF
    #define DISABLE_WARNING_USELESS_CAST
#endif

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define DO_PRAGMA(X)                              _Pragma(#X)                                    // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING_PUSH                      DO_PRAGMA(GCC diagnostic push)                 // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING_POP                       DO_PRAGMA(GCC diagnostic pop)                  // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING(warningName)              DO_PRAGMA(GCC diagnostic ignored warningName)  // NOLINT(cppcoreguidelines-macro-usage)

    #define DISABLE_WARNING_DEPRECATED_DECLARATIONS   DISABLE_WARNING("-Wdeprecated-declarations")   // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING_NULL_DEREF                DISABLE_WARNING("-Wnull-dereference")          // NOLINT(cppcoreguidelines-macro-usage)
#endif

// Defines for GCC only
#if defined(__GNUC__) && !defined(__clang__)
    #define DISABLE_WARNING_USELESS_CAST   DISABLE_WARNING("-Wuseless-cast")   // NOLINT(cppcoreguidelines-macro-usage)
#endif

// Defines for Clang only
#ifdef __clang__
    #define DISABLE_WARNING_USELESS_CAST
#endif

// Defines for unknown/unsupported compilers
#if !defined(_MSC_VER) && !defined(__GNUC__) && !defined(__clang__)
  #define DISABLE_WARNING
  #define DISABLE_WARNING_PUSH
  #define DISABLE_WARNING_POP

  #define DISABLE_WARNING_DEPRECATED_DECLARATIONS
  #define DISABLE_WARNING_NULL_DEREF
  #define DISABLE_WARNING_USELESS_CAST
#endif

// clang-format on
