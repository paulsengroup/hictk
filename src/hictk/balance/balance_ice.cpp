// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <filesystem>
#include <memory>

#include "./common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/balance.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

int run_subcmd(const BalanceICEConfig& c) {  // NOLINT(misc-use-internal-linkage)
  SPDLOG_INFO(FMT_STRING("balancing using ICE ({})"), c.name);
  const auto tmp_dir =
      !c.in_memory ? std::make_unique<const internal::TmpDir>(c.tmp_dir, true) : nullptr;
  const std::filesystem::path& tmp_dir_path = tmp_dir ? (*tmp_dir)() : "";

  if (hic::utils::is_hic_file(c.path_to_input.string())) {
    return balance_hic<balancing::ICE>(c, tmp_dir_path);
  }
  if (cooler::utils::is_multires_file(c.path_to_input.string())) {
    return balance_multires_cooler<balancing::ICE>(c, tmp_dir_path);
  }
  auto clr = cooler::File(c.path_to_input.string());
  return balance_cooler<balancing::ICE>(clr, c, tmp_dir_path);
}

}  // namespace hictk::tools
