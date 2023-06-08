// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// Source: https://www.fluentcpp.com/2019/08/30/how-to-disable-a-warning-in-cpp/

// clang-format off

// Defines for GCC and Clang
#if defined(__GNUC__) || defined(__clang__)
    #define DO_PRAGMA(X)                     _Pragma(#X)                                    // NOLINT(cppcoreguidelines-macro-usage)
    #define DISABLE_WARNING_PUSH             DO_PRAGMA(GCC diagnostic push)
    #define DISABLE_WARNING_POP              DO_PRAGMA(GCC diagnostic pop)
    #define DISABLE_WARNING(warningName)     DO_PRAGMA(GCC diagnostic ignored warningName)  // NOLINT(cppcoreguidelines-macro-usage)
#endif

// Defines specific to Clang
#ifdef __clang__
    #define DISABLE_WARNING_USELESS_CAST
#endif

// Defines specific to GCC
#if defined(__GNUC__) && !defined(__clang__)
    #define DISABLE_WARNING_USELESS_CAST  DISABLE_WARNING("-Wuseless-cast")
#endif

// Defines for unknown/unsupported compilers
#if !defined(__GNUC__) && !defined(__clang__)
    #define DISABLE_WARNING
    #define DISABLE_WARNING_PUSH
    #define DISABLE_WARNING_POP

    #define DISABLE_WARNING_USELESS_CAST
#endif

// clang-format on
