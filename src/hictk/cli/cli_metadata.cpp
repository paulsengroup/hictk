// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/validators.hpp"

namespace hictk::tools {

void Cli::make_metadata_subcommand() {
  auto& sc = *_cli.add_subcommand("metadata", "Print file metadata to stdout.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = MetadataConfig{};
                  });

  _config = MetadataConfig{};
  auto& c = std::get<MetadataConfig>(_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic or .[ms]cool file (Cooler URI syntax supported).")
      ->check(IsValidCoolerFile | IsValidHiCFile)
      ->required();

  sc.add_option(
      "-f,--output-format",
      c.output_format,
      "Format used to return file metadata.\n"
      "Should be one of: json, toml, or yaml.")
      ->check(CLI::IsMember({"json", "toml", "yaml"}))
      ->capture_default_str();
  sc.add_flag(
      "--include-file-path,!--exclude-file-path",
      c.include_file_path,
      "Output the given input path using attribute \"uri\".")
      ->capture_default_str();
  sc.add_flag(
      "--recursive",
      c.recursive,
      "Print metadata for each resolution or cell contained in a "
      "multi-resolution or single-cell file.")
      ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::transform_args_metadata_subcommand() {
  auto& c = std::get<MetadataConfig>(_config);

  c.input_format = infer_input_format(c.uri);

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
  c.verbosity = parse_hictk_verbosity_from_env().value_or(
      static_cast<std::int16_t>(spdlog::level::critical) - c.verbosity);
}

}  // namespace hictk::tools
