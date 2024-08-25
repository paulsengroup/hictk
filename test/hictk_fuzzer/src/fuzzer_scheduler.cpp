// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/fuzzer/config.hpp"
#include "hictk/fuzzer/tools.hpp"

namespace hictk::fuzzer {

int fuzz_subcommand([[maybe_unused]] const Config& c) { return 0; }

}  // namespace hictk::fuzzer
