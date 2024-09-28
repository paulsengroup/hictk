// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// NOLINTBEGIN(clang-analyzer-optin.core.EnumCastOutOfRange)
#if __has_include(<toml++/toml.hpp>)
#include <toml++/toml.hpp>
#elif __has_include(<toml++/toml.h>)
#include <toml++/toml.h>
#elif __has_include(<toml.hpp>)
#include <toml.hpp>
#else
#include <toml.h>
#endif
// NOLINTEND(clang-analyzer-optin.core.EnumCastOutOfRange)
