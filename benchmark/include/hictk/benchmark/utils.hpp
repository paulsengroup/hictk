// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <string>
#include <string_view>

#include "hictk/version.hpp"

namespace hictk::benchmark {

[[nodiscard]] inline std::string generate_test_name(std::string_view title,
                                                    bool add_braces = true) {
  return fmt::format(FMT_STRING("{}"
                                "\"name\": \"{}\", "
                                "\"hictk-version\": \"{}\", "
                                "\"start-time\": \"{:%FT%T}\""
                                "{}"),
                     add_braces ? "{" : "", title, config::version::str(),
                     fmt::gmtime(std::time(nullptr)), add_braces ? "}" : "");
}

}  // namespace hictk::benchmark
