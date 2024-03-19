// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

static void validate_base_resolution(const FixMcoolConfig& c, std::uint32_t base_resolution) {
  ValidateConfig vc{};
  vc.uri =
      fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_input.string(), base_resolution);
  vc.validate_index = true;

  [[maybe_unused]] const auto ec = validate_subcmd(vc);
  assert(ec == 0);
}

static void run_hictk_zoomify(const FixMcoolConfig& c,
                              const std::vector<std::uint32_t>& resolutions,
                              std::string_view base_uri) {
  ZoomifyConfig zc{};
  zc.path_to_input = std::string{base_uri};
  zc.path_to_output = c.path_to_output.string();
  zc.input_format = "cool";
  zc.output_format = "cool";
  zc.resolutions = resolutions;
  zc.copy_base_resolution = true;
  zc.force = c.force;
  zc.verbosity = c.verbosity;

  [[maybe_unused]] const auto ec = zoomify_subcmd(zc);
  assert(ec == 0);
}

static std::optional<BalanceICEConfig> detect_balancing_params(std::string_view file,
                                                               std::uint32_t resolution) {
  const HighFive::File clr(std::string{file}, HighFive::File::ReadOnly);
  const auto path = fmt::format(FMT_STRING("resolutions/{}/bins/weight"), resolution);

  if (!clr.exist(path)) {
    SPDLOG_WARN(
        FMT_STRING("Cooler at {}::{} does not appear to have been balanced. SKIPPING balancing!"),
        file, path);
    return {};
  }

  const cooler::Dataset dset(cooler::RootGroup{clr.getGroup("/")}, path);
  BalanceICEConfig c{};

  try {
    const auto cis_only = dset.read_attribute<bool>("cis_only");
    const auto trans_only =
        dset.has_attribute("trans_only") ? dset.read_attribute<bool>("trans_only") : false;

    assert(cis_only + trans_only < 2);

    if (cis_only) {
      c.mode = "cis";
      c.name = "ICE";
    } else if (trans_only) {
      c.mode = "trans";
      c.name = "INTER_ICE";
    } else {
      c.mode = "gw";
      c.name = "GW_ICE";
    }

    c.symlink_to_weight = true;

    c.masked_diags = dset.read_attribute<std::size_t>("ignore_diags");
    c.mad_max = dset.read_attribute<double>("mad_max");
    c.min_count = dset.read_attribute<std::size_t>("min_count");
    c.min_nnz = dset.read_attribute<std::size_t>("min_nnz");
    c.tolerance = dset.read_attribute<double>("tol");

    // NOLINTNEXTLINE
  } catch (const std::exception&) {
  }
  return c;
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

static void run_hictk_balance(const FixMcoolConfig& c, std::uint32_t resolution) {
  auto bc = detect_balancing_params(c.path_to_input.string(), resolution);

  if (!bc) {
    return;
  }

  bc->path_to_input = fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_output, resolution);
  bc->tmp_dir = c.tmp_dir;
  bc->in_memory = c.in_memory;
  bc->threads = c.threads;
  bc->zstd_compression_lvl = c.zstd_compression_lvl;
  bc->chunk_size = c.chunk_size;

  [[maybe_unused]] const auto ec = balance_subcmd(*bc);
  assert(ec == 0);
}

int fix_mcool_subcmd(const FixMcoolConfig& c) {
  assert(cooler::utils::is_multires_file(c.path_to_input.string()));

  const auto t0 = std::chrono::system_clock::now();

  const auto resolutions = cooler::MultiResFile{c.path_to_input.string()}.resolutions();
  const auto base_resolution = resolutions.front();

  const auto base_uri =
      fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_input.string(), base_resolution);

  if (c.check_base_resolution) {
    SPDLOG_INFO(FMT_STRING("Validating {}..."), base_uri);
    validate_base_resolution(c, base_resolution);
  }

  run_hictk_zoomify(c, resolutions, base_uri);

  std::for_each(resolutions.begin() + 1, resolutions.end(), [&](const auto& res) {
    run_hictk_balance(c, res);  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
  });

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

  SPDLOG_INFO(FMT_STRING("Restoration successfully completed! Elapsed time: {}s"),
              static_cast<double>(delta) / 1000.0);

  return 0;
}
}  // namespace hictk::tools
