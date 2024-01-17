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
#include <limits>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_merge_subcommand() {
  auto& sc =
      *_cli.add_subcommand("merge", "Merge multiple Cooler or .hic files into a single file.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = MergeConfig{};
           });

  _config = MergeConfig{};
  auto& c = std::get<MergeConfig>(_config);

  // clang-format off
  sc.add_option(
      "input-files",
      c.input_files,
      "Path to two or more Cooler or .hic files to be merged (Cooler URI syntax supported).")
      ->check(IsValidCoolerFile | IsValidHiCFile)
      ->expected(2, std::numeric_limits<int>::max())
      ->required();

  sc.add_option(
      "-o,--output-file",
      c.output_file,
      "Output Cooler or .hic file (Cooler URI syntax supported).")
      ->required();

  sc.add_option(
      "--resolution",
      c.resolution,
      "HiC matrix resolution (ignored when input files are in .cool format).")
      ->check(CLI::NonNegativeNumber);

  sc.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite output file.")
      ->capture_default_str();

  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Number of pixels to store in memory before writing to disk.")
      ->capture_default_str();

  sc.add_option(
      "-l,--compression-lvl",
      c.compression_lvl,
      "Compression level used to compress interactions.\n"
      "Defaults to 6 and 12 for .cool and .hic files, respectively.")
      ->check(CLI::Bound(1, 12));

  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.\n"
      "When merging interactions in Cooler format, only a single thread will be used.")
      ->check(CLI::Range(std::uint32_t(1), std::thread::hardware_concurrency()))
      ->capture_default_str();

  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->capture_default_str();

  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();

  // clang-format on

  _config = std::monostate{};
}

static bool check_all_files_are_hic(const std::vector<std::string>& paths) {
  return std::all_of(paths.begin(), paths.end(),
                     [](const auto& path) { return hic::utils::is_hic_file(path); });
}

static bool check_all_files_are_cooler(const std::vector<std::string>& paths) {
  return std::all_of(paths.begin(), paths.end(),
                     [](const auto& path) { return cooler::utils::is_cooler(path); });
}

static bool check_all_files_are_multires_cooler(const std::vector<std::string>& paths) {
  return std::all_of(paths.begin(), paths.end(),
                     [](const auto& path) { return cooler::utils::is_multires_file(path); });
}

static bool check_all_files_are_singlecell_cooler(const std::vector<std::string>& paths) {
  return std::all_of(paths.begin(), paths.end(),
                     [](const auto& path) { return cooler::utils::is_scool_file(path); });
}

void Cli::validate_merge_subcommand() const {
  assert(_cli.get_subcommand("merge")->parsed());

  std::vector<std::string> errors;
  std::vector<std::string> warnings;
  const auto& c = std::get<MergeConfig>(_config);
  const auto& sc = *_cli.get_subcommand("merge");

  if (!c.force && std::filesystem::exists(c.output_file)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_file));
  }

  const auto is_hic = check_all_files_are_hic(c.input_files);
  const auto is_cooler = check_all_files_are_cooler(c.input_files);
  const auto is_mcooler = check_all_files_are_multires_cooler(c.input_files);
  const auto is_scool = check_all_files_are_singlecell_cooler(c.input_files);

  if (is_scool) {
    errors.emplace_back("merging file in .scool format is not supported.");
  }

  const auto output_format = infer_output_format(c.output_file);

  auto input_output_format_mismatch = is_hic && output_format != "hic";
  input_output_format_mismatch |= ((is_cooler || is_mcooler) && output_format != "cool");

  if (input_output_format_mismatch) {
    errors.emplace_back(
        "detected mismatch in input-output formats: merging files of different formats is not "
        "supported.");
  }

  if (c.resolution == 0 && (is_hic || is_mcooler)) {
    errors.emplace_back("--resolution is mandatory when input files are in .hic or .mcool format.");
  }

  const auto resolution_parsed = !sc.get_option("--resolution")->empty();
  if (is_cooler && resolution_parsed) {
    warnings.emplace_back("--resolution is ignored when file is in .[s]cool format.");
  }

  for (const auto& w : warnings) {
    SPDLOG_WARN(FMT_STRING("{}"), w);
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_merge_subcommand() {
  auto& c = std::get<MergeConfig>(_config);
  const auto& sc = *_cli.get_subcommand("merge");

  c.output_format = c.output_file.empty() ? "text" : infer_output_format(c.output_file);

  if (sc.get_option("--compression-lvl")->empty()) {
    c.compression_lvl = c.output_format == "hic" ? 12 : 6;
  }

  c.tmp_dir /= (std::filesystem::path(c.output_file).filename().string() + ".tmp");

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
