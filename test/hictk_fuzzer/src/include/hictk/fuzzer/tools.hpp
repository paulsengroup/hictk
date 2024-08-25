// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/fuzzer/config.hpp"

namespace hictk::fuzzer {

[[nodiscard]] int fuzz_subcommand(const Config& c);
[[nodiscard]] int launch_worker_subcommand(const Config& c);

}  // namespace hictk::fuzzer
