// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <thread>
#include <variant>
#include <vector>

#include "hictk/file.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_balance_subcommand() {
  auto& sc = *_cli.add_subcommand("balance", "Balance HiC matrices using ICE or SCALE.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = BalanceConfig{};
                  });

  _config = BalanceConfig{};
  auto& c = std::get<BalanceConfig>(_config);

  // clang-format off
  sc.add_option(
      "input",
      c.path_to_input,
      "Path to the .hic, .cool or .mcool file to be balanced.")
      ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
      ->required();
  sc.add_flag(
     "--algorithm",
     c.algorithm,
     "Balancing algorithm.")
     ->check(CLI::IsMember({"ICE", "SCALE", "VC"}))
     ->capture_default_str();
  sc.add_option(
      "--mode",
      c.mode,
      "Balance matrix using:\n"
      " - genome-wide interactions (gw)\n"
      " - trans-only interactions (trans)\n"
      " - cis-only interactions (cis)")
      ->check(CLI::IsMember({"gw", "trans", "cis"}))
      ->capture_default_str();
  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->capture_default_str();
  sc.add_option(
      "--ignore-diags",
      c.masked_diags,
      "Number of diagonals (including the main diagonal) to mask before balancing.")
      ->capture_default_str();
  sc.add_option(
      "--mad-max",
      c.mad_max,
      "Mask bins using the MAD-max filter.\n"
      "bins whose log marginal sum is less than --mad-max median\n"
      "absolute deviations below the median log marginal sum of\n"
      "all the bins in the same chromosome.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_option(
      "--min-nnz",
      c.min_nnz,
      "Mask rows with fewer than --min-nnz non-zero entries.")
      ->capture_default_str();
  sc.add_option(
      "--min-count",
      c.min_count,
      "Mask rows with fewer than --min-count interactions.")
      ->capture_default_str();
  sc.add_option(
      "--tolerance",
      c.tolerance,
      "Threshold of the variance of marginals used to determine whether\n"
      "the algorithm has converged.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_option(
      "--max-iters",
      c.max_iters,
      "Maximum number of iterations.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();
  sc.add_flag(
      "--rescale-weights,!--no-rescale-weights",
      c.rescale_marginals,
      "Rescale weights such that rows sum approximately to 2.")
      ->capture_default_str();
  sc.add_option(
      "--name",
      c.name,
      "Name to use when writing weights to file.\n"
      "Defaults to ICE, INTER_ICE and GW_ICE when --mode is cis, trans and gw, respectively.")
      ->capture_default_str();
  sc.add_flag(
      "--create-weight-link" ,
      c.symlink_to_weight,
      "Create a symbolic link to the balancing weights at clr::/bins/weight.\n"
      "Ignored when balancing .hic files")
      ->capture_default_str();
  sc.add_flag(
      "--in-memory",
      c.in_memory,
      "Store all interactions in memory (greatly improves performance).")
      ->capture_default_str();
  sc.add_flag(
      "--stdout",
      c.stdout_,
      "Write balancing weights to stdout instead of writing them to the input file.")
      ->capture_default_str();
  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Number of interactions to process at once. Ignored when using --in-memory.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();
  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();
  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.")
      ->check(CLI::Range(std::uint32_t(1), std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_option(
      "-l,--compression-lvl",
      c.zstd_compression_lvl,
      "Compression level used to compress temporary files using ZSTD.")
      ->check(CLI::Range(0, 19))
      ->capture_default_str();
  sc.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files and datasets (if any).")
      ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::validate_balance_subcommand() const {
  [[maybe_unused]] const auto& c = std::get<BalanceConfig>(_config);
  std::vector<std::string> errors;

  const auto input_format = infer_input_format(c.path_to_input);
  if (input_format == "hic") {
    const auto avail_resolutions = hic::utils::list_resolutions(c.path_to_input);
    const hic::File f(c.path_to_input.string(), avail_resolutions.back());
    if (f.version() < 9) {
      errors.emplace_back("balancing .hic files v8 and older is not currently supported.");
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_balance_subcommand() {
  auto& c = std::get<BalanceConfig>(_config);

  if (c.name.empty()) {
    if (c.mode == "cis") {
      c.name = c.algorithm;
    } else if (c.mode == "trans") {
      c.name = "INTER_" + c.algorithm;
    } else {
      assert(c.mode == "gw");
      c.name = "GW_" + c.algorithm;
    }
  }

  const auto input_format = infer_input_format(c.path_to_input);
  auto input_path = c.path_to_input;
  if (input_format == "cool") {
    input_path = cooler::File(c.path_to_input.string()).path();
  }

  c.tmp_dir /= input_path.filename().string() + ".tmp";

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
