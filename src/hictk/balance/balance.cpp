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
#include "hictk/balancing/scale.hpp"
#include "hictk/balancing/vc.hpp"
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

template <typename BalanceConfig>
static void write_weights_hic(
    hic::internal::HiCFileWriter& hfw, const BalanceConfig& c,
    const phmap::flat_hash_map<std::uint32_t, balancing::Weights>& weights, bool force_overwrite) {
  for (const auto& [resolution, weights_] : weights) {
    hfw.add_norm_vector(c.name, "BP", resolution, weights_, force_overwrite);
  }
  hfw.write_norm_vectors_and_norm_expected_values();
}

template <typename BalanceConfig>
static void write_weights_cooler(std::string_view uri, const BalanceConfig& c,
                                 const balancing::Weights& weights,
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
  dset.append(weights.begin(), weights.end());

  dset.write_attribute("cis_only", c.mode == "cis");
  dset.write_attribute("divisive_weights", weights.type() == balancing::Weights::Type::DIVISIVE);
  if constexpr (std::is_same_v<BalanceICEConfig, BalanceConfig>) {
    dset.write_attribute("ignore_diags", std::int64_t(c.masked_diags));
    dset.write_attribute("mad_max", std::int64_t(c.mad_max));
    dset.write_attribute("min_count", std::int64_t(c.min_count));
    dset.write_attribute("min_nnz", std::int64_t(c.min_nnz));
    dset.write_attribute("tol", c.tolerance);

    if (c.mode != "cis") {
      if (variance.front() != -1) {
        dset.write_attribute("converged", variance.front() < c.tolerance);
        dset.write_attribute("scale", scale.front());
        dset.write_attribute("var", variance.front());
      }
    } else {
      std::vector<bool> converged{};
      for (const auto& var : variance) {
        if (var != -1) {
          converged.push_back(var < c.tolerance);  // NOLINT
        }
      }
      if (!converged.empty()) {
        dset.write_attribute("converged", converged);
        dset.write_attribute("scale", scale);
        dset.write_attribute("var", variance);
      }
    }
  }

  if (c.symlink_to_weight) {
    SPDLOG_INFO(FMT_STRING("Linking weights to {}::{}..."), file, link_path);
    if (clr.exist(link_path)) {
      clr.unlink(link_path);
    }
    clr.getGroup(grp).createSoftLink(link_path, dset());
  }
}

template <typename BalanceConfig>
static void write_weights_cooler(std::string_view uri, const BalanceConfig& c,
                                 const balancing::Weights& weights) {
  return write_weights_cooler(uri, c, weights, {-1}, {-1});
}

template <typename Balancer, typename BalanceConfig>
static auto init_params(const BalanceConfig& c, const std::filesystem::path& tmpfile) ->
    typename Balancer::Params {
  if constexpr (std::is_same_v<Balancer, balancing::ICE>) {
    return {c.tolerance, c.max_iters, c.masked_diags, c.min_nnz, c.min_count,
            c.mad_max,   tmpfile,     c.chunk_size,   c.threads};
  }

  if constexpr (std::is_same_v<Balancer, balancing::SCALE>) {
    return {c.tolerance, c.max_iters, 10.0, 1.0e-5, 0.05, 0.05, tmpfile, c.chunk_size, c.threads};
  }

  return {};
}

template <typename Balancer, typename BalanceConfig>
static int balance_cooler(cooler::File& f, const BalanceConfig& c,
                          const std::filesystem::path& tmp_dir) {
  if (!c.force && !c.stdout_ && f.has_normalization(c.name)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Normalization weights for \"{}\" already exist in file {}. Pass "
                               "--force to overwrite existing weights."),
                    c.name, f.path()));
  }

  const auto tmpfile =
      tmp_dir.empty() ? ""
                      : tmp_dir / (std::filesystem::path{f.path()}.filename().string() + ".tmp");
  const auto params = init_params<Balancer>(c, tmpfile);

  typename Balancer::Type mode{};
  if (c.mode == "gw") {
    mode = Balancer::Type::gw;
  } else if (c.mode == "cis") {
    mode = Balancer::Type::cis;
  } else {
    mode = Balancer::Type::trans;
  }

  const Balancer balancer(f, mode, params);
  const auto weights = balancer.get_weights(c.rescale_marginals);

  if (c.stdout_) {
    for (const auto& w :
         balancer.get_weights(c.rescale_marginals)(balancing::Weights::Type::DIVISIVE)) {
      fmt::print(FMT_COMPILE("{}\n"), w);
    }
    return 0;
  }

  const auto uri = f.uri();
  f.close();
  if constexpr (std::is_same_v<Balancer, balancing::ICE>) {
    write_weights_cooler(uri, c, weights, balancer.variance(), balancer.scale());
  } else {
    write_weights_cooler(uri, c, weights);
  }
  return 0;
}

template <typename Balancer, typename BalanceConfig>
static int balance_hic(const BalanceConfig& c, const std::filesystem::path& tmp_dir) {
  const auto resolutions = hic::utils::list_resolutions(c.path_to_input);
  for (const auto& res : resolutions) {
    const hic::File f(c.path_to_input.string(), res);
    if (!c.force && !c.stdout_ && f.has_normalization(c.name)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Normalization weights for \"{}\" already exist in file {}. Pass "
                                 "--force to overwrite existing weights."),
                      c.name, f.path()));
    }
  }

  const auto tmpfile =
      tmp_dir / (std::filesystem::path{c.path_to_input}.filename().string() + ".tmp");

  const auto params = init_params<Balancer>(c, tmpfile);
  typename Balancer::Type mode{};
  if (c.mode == "gw") {
    mode = Balancer::Type::gw;
  } else if (c.mode == "cis") {
    mode = Balancer::Type::cis;
  } else {
    mode = Balancer::Type::trans;
  }

  phmap::flat_hash_map<std::uint32_t, balancing::Weights> weights{resolutions.size()};
  for (const auto& res : resolutions) {
    SPDLOG_INFO(FMT_STRING("balancing resolution {}..."), res);
    const hic::File f(c.path_to_input.string(), res);

    const Balancer balancer(f, mode, params);

    if (c.stdout_) {
      for (const auto& w :
           balancer.get_weights(c.rescale_marginals)(balancing::Weights::Type::DIVISIVE)) {
        fmt::print(FMT_COMPILE("{}\n"), w);
      }
      return 0;
    }
    weights.emplace(res, balancer.get_weights(c.rescale_marginals));
  }

  // NOLINTNEXTLINE(misc-const-correctness)
  hic::internal::HiCFileWriter hfw(c.path_to_input.string());
  write_weights_hic(hfw, c, weights, c.force);
  return 0;
}

template <typename Balancer, typename BalanceConfig>
static int balance_multires_cooler(const BalanceConfig& c, const std::filesystem::path& tmp_dir) {
  const auto resolutions = cooler::utils::list_resolutions(c.path_to_input.string());

  for (const auto& res : resolutions) {
    auto clr = cooler::MultiResFile(c.path_to_input.string()).open(res);
    SPDLOG_INFO(FMT_STRING("balancing resolution {}..."), res);
    balance_cooler<Balancer>(clr, c, tmp_dir);
  }
  return 0;
}

int balance_subcmd(const BalanceICEConfig& c) {
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

int balance_subcmd(const BalanceSCALEConfig& c) {
  SPDLOG_INFO(FMT_STRING("balancing using SCALE ({})"), c.name);
  const auto tmp_dir =
      !c.in_memory ? std::make_unique<const internal::TmpDir>(c.tmp_dir, true) : nullptr;
  const std::filesystem::path& tmp_dir_path = tmp_dir ? (*tmp_dir)() : "";

  if (hic::utils::is_hic_file(c.path_to_input.string())) {
    return balance_hic<balancing::SCALE>(c, tmp_dir_path);
  }
  if (cooler::utils::is_multires_file(c.path_to_input.string())) {
    return balance_multires_cooler<balancing::SCALE>(c, tmp_dir_path);
  }
  auto clr = cooler::File(c.path_to_input.string());
  return balance_cooler<balancing::SCALE>(clr, c, tmp_dir_path);
}

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
