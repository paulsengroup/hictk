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
  auto& sc = *_cli.add_subcommand("merge", "Merge coolers.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = MergeConfig{};
                  });

  _config = MergeConfig{};
  auto& c = std::get<MergeConfig>(_config);

  // clang-format off
  sc.add_option(
      "input-coolers",
      c.input_uris,
      "Path to two or more Cooler files to be merged (URI syntax supported).")
      ->check(IsValidCoolerFile)
      ->expected(2, std::numeric_limits<int>::max())
      ->required();

  sc.add_option(
      "-o,--output-cooler",
      c.output_uri,
      "Output Cooler (URI syntax supported).\n"
      "When not specified, merged interactions will be printed to stdout.");

  sc.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite output cooler.")
      ->capture_default_str();

  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Number of pixels to store in memory before writing to disk.")
      ->capture_default_str();

  sc.add_option(
      "-l,--compression-lvl",
      c.compression_lvl,
      "Compression level used to compress interactions.")
      ->check(CLI::Bound(1, 9))
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

void Cli::validate_merge_subcommand() const {
  assert(_cli.get_subcommand("merge")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<MergeConfig>(_config);

  if (!c.force && std::filesystem::exists(c.output_uri)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_uri));
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

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
