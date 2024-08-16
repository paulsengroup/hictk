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
#include <variant>
#include <vector>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_load_subcommand() {
  auto& sc =
      *_cli.add_subcommand("load",
                           "Build .cool and .hic files from interactions in various text formats.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = LoadConfig{};
           });

  _config = LoadConfig{};
  auto& c = std::get<LoadConfig>(_config);

  // clang-format off
  sc.add_option(
      "chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to .chrom.sizes file.")
      ->check(CLI::ExistingFile)
      ->required();

  sc.add_option(
      "output-path",
      c.output_path,
      "Path to output file.")
      ->required();

  sc.add_option(
      "-b,--bin-size",
      c.bin_size,
      "Bin size (bp).\n"
      "Required when --bin-table is not used.")
      ->check(CLI::PositiveNumber);

  sc.add_option(
      "--bin-table",
      c.path_to_bin_table,
      "Path to a BED3+ file with the bin table.")
      ->check(CLI::ExistingFile);

  sc.add_option(
      "-f,--format",
      c.format,
      "Input format.")
      ->check(CLI::IsMember({"4dn", "validpairs", "bg2", "coo"}))
      ->required();

  sc.add_flag(
      "--force",
      c.force,
      "Force overwrite existing output file(s).")
      ->capture_default_str();

  sc.add_option(
      "--assembly",
      c.assembly,
      "Assembly name.")
      ->capture_default_str();

  sc.add_flag(
      "--one-based,!--zero-based",
      c.one_based,
      "Interpret genomic coordinates or bins as one/zero based.\n"
      "By default coordinates are assumed to be one-based for interactions in\n"
      "4dn and validapairs formats and zero-based otherwise.");

  sc.add_flag(
      "--count-as-float",
      c.count_as_float,
      "Interactions are floats.")
      ->capture_default_str();

  sc.add_flag(
      "--skip-all-vs-all,!--no-skip-all-vs-all",
      c.skip_all_vs_all_matrix,
      "Do not generate All vs All matrix.\n"
      "Has no effect when creating .cool files.")
      ->capture_default_str();

  sc.add_flag(
      "--assume-sorted,!--assume-unsorted",
      c.assume_sorted,
      "Assume input files are already sorted.")
      ->capture_default_str();

  sc.add_option(
      "--chunk-size",
      c.batch_size,
      "Number of pixels to buffer in memory.")
      ->capture_default_str();

  sc.add_option(
      "-l,--compression-lvl",
      c.compression_lvl,
      "Compression level used to compress interactions.\n"
      "Defaults to 6 and 10 for .cool and .hic files, respectively.")
      ->check(CLI::Bound(1, 12));

  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.\n"
      "When loading interactions in a .cool file, only a single thread will be used.")
      ->check(CLI::Range(std::uint32_t(1), std::thread::hardware_concurrency()))
      ->capture_default_str();

  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->check(CLI::ExistingDirectory)
      ->capture_default_str();


  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();
  // clang-format on

  sc.get_option("--bin-size")->excludes(sc.get_option("--bin-table"));
  _config = std::monostate{};
}

void Cli::validate_load_subcommand() const {
  assert(_cli.get_subcommand("load")->parsed());

  std::vector<std::string> warnings;
  std::vector<std::string> errors;
  const auto& c = std::get<LoadConfig>(_config);
  const auto& sc = *_cli.get_subcommand("load");

  if (!c.force && std::filesystem::exists(c.output_path)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_path));
  }

  if (c.path_to_bin_table.empty() && c.bin_size == 0) {
    assert(c.bin_size == 0);
    errors.emplace_back("--bin-size is required when --bin-table is not specified.");
  }

  const auto output_format = infer_output_format(c.output_path);
  if (!c.path_to_bin_table.empty() && output_format == "hic") {
    errors.emplace_back("--bin-table is not supported when generating .hic files.");
  }

  if ((c.format == "bg2" || c.format == "coo") && !sc.get_option("--bin-table")->empty()) {
    errors.emplace_back(
        "specifying bins through the --bin-table is not supported when ingesting pre-binned "
        "interactions.");
  }

  if (c.format == "4dn" && c.format == "validpairs" && c.assume_sorted) {
    warnings.emplace_back(
        "--assume-sorted has no effect when ingesting interactions in 4dn or validpairs format.");
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

void Cli::transform_args_load_subcommand() {
  auto& c = std::get<LoadConfig>(_config);
  const auto& sc = *_cli.get_subcommand("load");

  c.output_format = infer_output_format(c.output_path);

  if (sc.get_option("--one-based")->empty()) {
    if (c.format == "4dn" || c.format == "validpairs") {
      c.offset = -1;
    }
  } else {
    c.offset = c.one_based ? -1 : 0;
  }

  if (sc.get_option("--compression-lvl")->empty()) {
    c.compression_lvl = c.output_format == "hic" ? 10 : 6;
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
