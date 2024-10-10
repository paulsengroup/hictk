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

#include "hictk/file.hpp"
#include "hictk/tmpdir.hpp"
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
      "--output-fmt",
      c.output_format,
      "Output format (by default this is inferred from the output file extension).\n"
      "Should be one of:\n"
      "- cool\n"
      "- hic\n")
      ->check(CLI::IsMember({"cool", "hic"}))
      ->default_str("auto");
  sc.add_option(
      "--resolution",
      c.resolution,
      "Hi-C matrix resolution (ignored when input files are in .cool format).")
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
      "Defaults to 6 and 10 for .cool and .hic files, respectively.")
      ->check(CLI::Bound(std::uint8_t{1}, MAX_HIC_COMPRESSION_LEVEL));
  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.\n"
      "When merging interactions in Cooler format, only a single thread will be used.")
      ->check(CLI::Range(std::uint32_t{1}, std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->check(CLI::ExistingDirectory)
      ->capture_default_str();
  sc.add_flag(
      "--skip-all-vs-all,!--no-skip-all-vs-all",
      c.skip_all_vs_all_matrix,
      "Do not generate All vs All matrix.\n"
      "Has no effect when merging .cool files.")
      ->capture_default_str();
  sc.add_option(
      "--count-type",
      c.count_type,
      "Specify the count type to be used when merging files.\n"
      "Ignored when the output file is in .hic format.")
      ->check(CLI::IsMember{{"int", "float"}})
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

static void validate_files_format(const std::vector<std::string>& paths, std::uint32_t resolution,
                                  std::vector<std::string>& errors) {
  assert(!paths.empty());
  const auto& p0 = paths.front();
  if ((hic::utils::is_hic_file(p0) || cooler::utils::is_multires_file(p0)) && resolution == 0) {
    errors.emplace_back("--resolution is mandatory when input files are in .hic or .mcool format.");
    return;
  }

  for (const auto& p : paths) {
    const auto format = infer_input_format(p);
    if (format == "scool") {
      errors.emplace_back("merging file in .scool format is not supported.");
      return;
    }
    if ((format == "hic" || format == "mcool") && resolution == 0) {
      errors.emplace_back(
          "--resolution is mandatory when one or more input files are in .hic or .mcool format.");
      return;
    }
  }
}

void Cli::validate_merge_subcommand() const {
  assert(_cli.get_subcommand("merge")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<MergeConfig>(_config);

  if (!c.force && std::filesystem::exists(c.output_file)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_file));
  }

  validate_files_format(c.input_files, c.resolution, errors);

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

  c.output_format = infer_output_format(c.output_file);

  if (c.resolution == 0) {
    c.resolution = File(c.input_files.front()).resolution();
  }

  if (sc.get_option("--compression-lvl")->empty()) {
    c.compression_lvl =
        c.output_format == "hic" ? DEFAULT_HIC_COMPRESSION_LEVEL : DEFAULT_COOL_COMPRESSION_LEVEL;
  }

  if (sc.get_option("--tmpdir")->empty()) {
    c.tmp_dir = internal::TmpDir::default_temp_directory_path();
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
