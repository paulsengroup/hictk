// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <spdlog/spdlog.h>

#include <variant>

#include "hictk/balancing/ice.hpp"
#include "hictk/balancing/methods.hpp"
#include "hictk/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/tools/common.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/juicer_tools.hpp"

namespace hictk::tools {

static void write_weights_hic(const hic::File& hf, const BalanceConfig& c,
                              const std::vector<double>& weights) {
  auto tmpfile = c.tmp_dir / std::filesystem::path{hf.name()}.filename();
  for (std::size_t i = 0; i < 1024; ++i) {
    if (!std::filesystem::exists(tmpfile)) {
      break;
    }

    tmpfile.replace_extension(".tmp" + std::to_string(i));
  }

  if (std::filesystem::exists(tmpfile)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to create temporary file {}"), tmpfile));
  }

  try {
    {
      const std::unique_ptr<FILE> f(std::fopen(tmpfile.string().c_str(), "ae"));
      if (!bool(f)) {
        throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), tmpfile);
      }

      std::ptrdiff_t i0 = 0;

      for (const auto& chrom : hf.chromosomes()) {
        if (chrom.is_all()) {
          continue;
        }
        fmt::print(f.get(), FMT_STRING("vector\t{}\t{}\t{}\tBP\n"), c.name, chrom.name(),
                   hf.bin_size());

        const auto num_bins = (chrom.size() + hf.bin_size() - 1) / hf.bin_size();
        const auto i1 = i0 + static_cast<std::ptrdiff_t>(num_bins);
        std::for_each(weights.begin() + i0, weights.begin() + i1, [&](const double w) {
          std::isnan(w) ? fmt::print(f.get(), FMT_COMPILE(".\n"))
                        : fmt::print(f.get(), FMT_COMPILE("{}\n"), 1.0 / w);
          if (!bool(f)) {  // NOLINT
            throw fmt::system_error(
                errno, FMT_STRING("an error occurred while writing weights to file {}"), tmpfile);
          }
        });

        i0 = i1;
      }
    }

    auto jt = run_juicer_tools_add_norm(c.juicer_tools_jar, tmpfile, hf.url(), c.juicer_tools_xmx);
    jt->wait();
    if (jt->exit_code() != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("juicer_tools pre failed with exit code {}"), jt->exit_code()));
    }
  } catch (...) {
    std::error_code ec{};
    std::filesystem::remove(tmpfile, ec);
  }
  std::filesystem::remove(tmpfile);
}

static void write_weights_cooler(std::string_view uri, const BalanceConfig& c,
                                 const std::vector<double>& weights,
                                 const std::vector<double>& variance,
                                 const std::vector<double>& scale) {
  const auto& [file, grp] = cooler::parse_cooler_uri(uri);
  const auto path = fmt::format(FMT_STRING("{}/bins/{}"), grp, c.name);
  SPDLOG_INFO(FMT_STRING("Writing weights to {}{}..."), uri, path);

  const HighFive::File clr(file, HighFive::File::ReadWrite);

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
      converged.push_back(var < c.tolerance);
    }
    dset.write_attribute("converged", converged);
    dset.write_attribute("scale", scale);
    dset.write_attribute("var", variance);
  }
}

static int balance_singleres_file(File&& f, const BalanceConfig& c) {
  std::filesystem::path tmpfile{};

  if (f.is_cooler()) {
    const auto& ff = f.get<cooler::File>();
    if (ff.has_weights(c.name) && !c.force) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "{}/bins/weight already exists. Pass --force to overwrite currently stored weights."),
          ff.uri()));
    }
  }

  if (!c.in_memory) {
    tmpfile = c.tmp_dir / std::filesystem::path{f.path()}.filename();
    for (std::size_t i = 0; i < 1024; ++i) {
      if (!std::filesystem::exists(tmpfile)) {
        break;
      }

      tmpfile.replace_extension(".tmp" + std::to_string(i));
    }

    if (std::filesystem::exists(tmpfile)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to create temporary file {}"), tmpfile));
    }
  }

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

  const auto balancer =
      std::visit([&](const auto& ff) { return balancing::ICE(ff, mode, params); }, f.get());
  const auto weights = balancer.get_weights(c.rescale_marginals);

  if (c.stdout_) {
    std::for_each(weights.begin(), weights.end(),
                  [&](const auto w) { fmt::print(FMT_COMPILE("{}\n"), w); });
    return 0;
  }

  if (f.is_cooler()) {
    const auto uri = f.uri();
    f.get<cooler::File>().close();
    write_weights_cooler(uri, c, weights, balancer.variance(), balancer.scale());
    return 0;
  }

  write_weights_hic(f.get<hic::File>(), c, weights);
  // TODO write weights .hic

  return 0;
}

static int balance_multires(const BalanceConfig& c) {
  const auto resolutions = cooler::MultiResFile(c.path_to_input.string()).resolutions();

  for (const auto& res : resolutions) {
    balance_singleres_file(
        File(fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_input.string(), res)), c);
  }
  return 0;
}

int balance_subcmd(const BalanceConfig& c) {
  if (cooler::utils::is_multires_file(c.path_to_input.string())) {
    return balance_multires(c);
  }

  std::vector<std::uint32_t> resolutions{};
  if (hic::utils::is_hic_file(c.path_to_input)) {
    resolutions = hic::utils::list_resolutions(c.path_to_input);
  } else {
    resolutions.push_back(File(c.path_to_input.string()).bin_size());
  }

  for (const auto& res : resolutions) {
    balance_singleres_file(File(c.path_to_input, res), c);
  }

  return 0;
}
}  // namespace hictk::tools
