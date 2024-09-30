// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "./common.hpp"
#include "hictk/balancing/vc.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

int balance_subcmd(const BalanceVCConfig& c) {
  SPDLOG_INFO(FMT_STRING("balancing using VC ({})"), c.name);
  if (hic::utils::is_hic_file(c.path_to_input.string())) {
    return balance_hic<balancing::VC>(c, "");
  }
  if (cooler::utils::is_multires_file(c.path_to_input.string())) {
    return balance_multires_cooler<balancing::VC>(c, "");
  }
  auto clr = cooler::File(c.path_to_input.string());
  return balance_cooler<balancing::VC>(clr, c, "");
}

}  // namespace hictk::tools
