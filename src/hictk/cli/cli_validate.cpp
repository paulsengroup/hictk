// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <cstddef>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {
void Cli::make_validate_subcommand() {
  auto& sc = *_cli.add_subcommand("validate", "Validate .hic and Cooler files.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = ValidateConfig{};
                  });

  _config = ValidateConfig{};
  auto& c = std::get<ValidateConfig>(_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic or .[ms]cool file (Cooler URI syntax supported).")
      ->required();
  sc.add_flag(
      "--validate-index",
      c.validate_index,
      "Validate Cooler index (may take a long time).")
      ->capture_default_str();
  sc.add_flag(
      "--validate-pixels",
      c.validate_pixels,
      "Validate pixels found in Cooler files (may take a long time).")
      ->capture_default_str();
  sc.add_option(
      "-f,--output-format",
      c.output_format,
      "Format used to report the outcome of file validation.\n"
      "Should be one of: json, toml, or yaml.")
      ->check(CLI::IsMember({"json", "toml", "yaml"}))
      ->capture_default_str();
  sc.add_flag(
      "--include-file-path,!--exclude-file-path",
      c.include_file_path,
      "Output the given input path using attribute \"uri\".")
      ->capture_default_str();
  sc.add_flag(
      "--exhaustive,!--fail-fast",
      c.exhaustive,
      "When processing multi-resolution or single-cell files,\n"
      "do not fail as soon as the first error is detected.")
      ->capture_default_str();
  sc.add_flag(
      "--quiet",
      c.quiet,
      "Don't print anything to stdout. Success/failure is reported through exit codes")
      ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::transform_args_validate_subcommand() {
  assert(_cli.get_subcommand("validate")->parsed());
  auto& c = std::get<ValidateConfig>(_config);

  if (c.quiet) {
    c.verbosity = static_cast<std::int16_t>(spdlog::level::err);
  } else {
    // in spdlog, high numbers correspond to low log levels
    assert(c.verbosity > 0 && c.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
    c.verbosity = parse_hictk_verbosity_from_env().value_or(
        static_cast<std::int16_t>(spdlog::level::critical) - c.verbosity);
  }
}

}  // namespace hictk::tools
