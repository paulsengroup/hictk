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
#include <filesystem>
#include <stdexcept>
#include <string>
#include <thread>
#include <variant>
#include <vector>

#include "hictk/tmpdir.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_fix_mcool_subcommand() {
  auto& sc = *_cli.add_subcommand("fix-mcool", "Fix corrupted .mcool files.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = FixMcoolConfig{};
                  });

  _config = FixMcoolConfig{};
  auto& c = std::get<FixMcoolConfig>(_config);

  // clang-format off
  sc.add_option(
      "input",
      c.path_to_input,
      "Path to a corrupted .mcool file.")
      ->check(IsValidMultiresCoolerFile)
      ->required();
  sc.add_option(
      "output",
      c.path_to_output,
      "Path where to store the restored .mcool.")
      ->required();
  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->check(CLI::ExistingDirectory)
      ->capture_default_str();
  sc.add_flag(
      "--skip-balancing",
      c.skip_balancing,
      "Do not recompute or copy balancing weights.");
  sc.add_flag(
      "--check-base-resolution",
      c.check_base_resolution,
      "Check whether the base resolution is corrupted.");
  sc.add_flag(
      "--in-memory",
      c.in_memory,
      "Store all interactions in memory while balancing (greatly improves performance).")
      ->capture_default_str();
  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Number of interactions to process at once during balancing.\n"
      "Ignored when using --in-memory.")
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
      "Maximum number of parallel threads to spawn (only applies to the balancing stage).")
      ->check(CLI::Range(std::uint32_t(1), std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_option(
      "-l,--compression-lvl",
      c.zstd_compression_lvl,
      "Compression level used to compress temporary files using ZSTD (only applies to the balancing stage).")
      ->check(CLI::Range(0, 19))
      ->capture_default_str();
  sc.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files (if any).")
      ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::validate_fix_mcool_subcommand() const {
  const auto& c = std::get<FixMcoolConfig>(_config);
  std::vector<std::string> errors;
  std::vector<std::string> warnings{};

  if (!c.force && std::filesystem::exists(c.path_to_output)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.path_to_output));
  }

  if (c.skip_balancing) {
    const auto* sc = _cli.get_subcommand("fix-mcool");
    if (!sc->get_option("--tmpdir")->empty()) {
      warnings.emplace_back("option --tmpdir is ignored when --skip-balancing is provided.");
    }
    if (!sc->get_option("--in-memory")->empty()) {
      warnings.emplace_back("option --in-memory is ignored when --skip-balancing is provided.");
    }
    if (!sc->get_option("--compression-lvl")->empty()) {
      warnings.emplace_back(
          "option --compression-lvl is ignored when --skip-balancing is provided.");
    }
    if (!sc->get_option("--chunk-size")->empty()) {
      warnings.emplace_back("option --chunk-size is ignored when --skip-balancing is provided.");
    }
    if (!sc->get_option("--threads")->empty()) {
      warnings.emplace_back("option --threads is ignored when --skip-balancing is provided.");
    }
  }

  for (const auto& w : warnings) {
    SPDLOG_WARN(FMT_STRING("{}"), w);
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_fix_mcool_subcommand() {
  auto& c = std::get<FixMcoolConfig>(_config);
  const auto& sc = *_cli.get_subcommand("fix-mcool");

  if (sc.get_option("--tmpdir")->empty()) {
    c.tmp_dir = hictk::internal::TmpDir::default_temp_directory_path();
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
