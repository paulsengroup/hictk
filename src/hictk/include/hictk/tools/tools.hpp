// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <stdexcept>
#include <variant>

#include "hictk/tools/balance.hpp"
#include "hictk/tools/convert.hpp"
#include "hictk/tools/dump.hpp"
#include "hictk/tools/fix_mcool.hpp"
#include "hictk/tools/load.hpp"
#include "hictk/tools/merge.hpp"
#include "hictk/tools/metadata.hpp"
#include "hictk/tools/rename_chromosomes.hpp"
#include "hictk/tools/validate.hpp"
#include "hictk/tools/zoomify.hpp"

namespace hictk::tools {

[[noreturn]] inline int run_subcmd([[maybe_unused]] const std::monostate& c) {
  throw std::logic_error(
      "hictk::tools::run_subcmd(const std::monostate& c) should never be called!");
}

}  // namespace hictk::tools
