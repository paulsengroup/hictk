// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/balancing/ice.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void write_weights_hic(
    hic::internal::HiCFileWriter& hfw, const BalanceConfig& c,
    const phmap::flat_hash_map<std::uint32_t, std::vector<double>>& weights, bool force_overwrite) {
  for (const auto& [resolution, weights_] : weights) {
    std::vector<float> weights_f(weights_.size());
    std::transform(weights_.begin(), weights_.end(), weights_f.begin(),
                   [](const auto w) { return static_cast<float>(1.0 / w); });
    hfw.add_norm_vector(c.name, "BP", resolution, weights_f, force_overwrite);
  }
  hfw.write_norm_vectors_and_norm_expected_values();
}

static void write_weights_cooler(std::string_view uri, const BalanceConfig& c,
                                 const std::vector<double>& weights,
                                 const std::vector<double>& variance,
                                 const std::vector<double>& scale) {
  const auto& [file, grp] = cooler::parse_cooler_uri(uri);
  const auto path = fmt::format(FMT_STRING("{}/bins/{}"), grp, c.name);
  const auto link_path = fmt::format(FMT_STRING("{}/bins/weight"), grp);

  SPDLOG_INFO(FMT_STRING("Writing weights to {}::{}..."), file, path);
  const HighFive::File clr(file, HighFive::File::ReadWrite);

  if (c.symlink_to_weight && clr.exist(link_path) && !c.force) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to create link to {}::{}: object already exists"), file, link_path));
  }

  if (clr.exist(path)) {
    assert(c.force);
    clr.unlink(path);
  }

  cooler::Dataset dset(cooler::RootGroup{clr.getGroup(grp)}, path, 0.0);
  dset.append(weights);

  dset.write_attribute("cis_only", c.mode == "cis");
  dset.write_attribute("divisive_weights", false);
  dset.write_attribute("ignore_diags", std::int64_t(c.masked_diags));
  dset.write_attribute("mad_max", std::int64_t(c.mad_max));
  dset.write_attribute("min_count", std::int64_t(c.min_count));
  dset.write_attribute("min_nnz", std::int64_t(c.min_nnz));
  dset.write_attribute("tol", c.tolerance);

  if (c.mode != "cis") {
    dset.write_attribute("converged", variance.front() < c.tolerance);
    dset.write_attribute("scale", scale.front());
    dset.write_attribute("var", variance.front());
  } else {
    std::vector<bool> converged{};
    for (const auto& var : variance) {
      converged.push_back(var < c.tolerance);  // NOLINT
    }
    dset.write_attribute("converged", converged);
    dset.write_attribute("scale", scale);
    dset.write_attribute("var", variance);
  }

  if (c.symlink_to_weight) {
    SPDLOG_INFO(FMT_STRING("Linking weights to {}::{}..."), file, link_path);
    if (clr.exist(link_path)) {
      clr.unlink(link_path);
    }
    clr.getGroup(grp).createSoftLink(link_path, dset());
  }
}

// NOLINTNEXTLINE(*-rvalue-reference-param-not-moved)
static int balance_cooler(cooler::File&& f, const BalanceConfig& c) {
  if (!c.force && !c.stdout_ && f.has_normalization(c.name)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Normalization weights for \"{}\" already exist in file {}. Pass "
                               "--force to overwrite existing weights."),
                    c.name, f.path()));
  }

  const auto tmpfile = c.tmp_dir / std::filesystem::path{f.path()}.filename();

  const balancing::ICE::Params params{c.tolerance, c.max_iters,  c.masked_diags,
                                      c.min_nnz,   c.min_count,  c.mad_max,
                                      tmpfile,     c.chunk_size, c.threads};
  balancing::ICE::Type mode{};
  if (c.mode == "gw") {
    mode = balancing::ICE::Type::gw;
  } else if (c.mode == "cis") {
    mode = balancing::ICE::Type::cis;
  } else {
    mode = balancing::ICE::Type::trans;
  }

  const balancing::ICE balancer(f, mode, params);
  const auto weights = balancer.get_weights(c.rescale_marginals);

  if (c.stdout_) {
    std::for_each(weights.begin(), weights.end(),
                  [&](const auto w) { fmt::print(FMT_COMPILE("{}\n"), w); });
    return 0;
  }

  const auto uri = f.uri();
  f.close();
  write_weights_cooler(uri, c, weights, balancer.variance(), balancer.scale());
  return 0;
}

// NOLINTNEXTLINE(*-rvalue-reference-param-not-moved)
static int balance_hic(const BalanceConfig& c) {
  const auto resolutions = hic::utils::list_resolutions(c.path_to_input);
  for (const auto& res : resolutions) {
    const hic::File f(c.path_to_input, res);
    if (!c.force && !c.stdout_ && f.has_normalization(c.name)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Normalization weights for \"{}\" already exist in file {}. Pass "
                                 "--force to overwrite existing weights."),
                      c.name, f.path()));
    }
  }

  const auto tmpfile = c.tmp_dir / std::filesystem::path{c.path_to_input}.filename();

  const balancing::ICE::Params params{c.tolerance, c.max_iters,  c.masked_diags,
                                      c.min_nnz,   c.min_count,  c.mad_max,
                                      tmpfile,     c.chunk_size, c.threads};
  balancing::ICE::Type mode{};
  if (c.mode == "gw") {
    mode = balancing::ICE::Type::gw;
  } else if (c.mode == "cis") {
    mode = balancing::ICE::Type::cis;
  } else {
    mode = balancing::ICE::Type::trans;
  }

  phmap::flat_hash_map<std::uint32_t, std::vector<double>> weights{resolutions.size()};
  for (const auto& res : resolutions) {
    SPDLOG_INFO(FMT_STRING("balancing resolution {}..."), res);
    const hic::File f(c.path_to_input, res);
    const balancing::ICE balancer(f, mode, params);

    if (c.stdout_) {
      std::for_each(weights.begin(), weights.end(),
                    [&](const auto w) { fmt::print(FMT_COMPILE("{}\n"), w); });
    } else {
      weights.emplace(res, balancer.get_weights(c.rescale_marginals));
    }
  }

  hic::internal::HiCFileWriter hfw(c.path_to_input.string(), c.threads);
  write_weights_hic(hfw, c, weights, c.force);
  return 0;
}

static int balance_multires_cooler(const BalanceConfig& c) {
  const cooler::MultiResFile mclr(c.path_to_input.string());

  for (const auto& res : mclr.resolutions()) {
    SPDLOG_INFO(FMT_STRING("balancing resolution {}..."), res);
    balance_cooler(mclr.open(res), c);
  }
  return 0;
}

int balance_subcmd(const BalanceConfig& c) {
  [[maybe_unused]] const internal::TmpDir tmp_dir{c.tmp_dir};

  if (hic::utils::is_hic_file(c.path_to_input.string())) {
    return balance_hic(c);
  }

  if (cooler::utils::is_multires_file(c.path_to_input.string())) {
    return balance_multires_cooler(c);
  }

  balance_cooler(cooler::File(c.path_to_input.string()), c);

  return 0;
}
}  // namespace hictk::tools
