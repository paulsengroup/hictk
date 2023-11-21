// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <string>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {
void Cli::make_rename_chromosomes_subcommand() {
  auto& sc =
      *_cli.add_subcommand("rename-chromosomes", "Rename chromosomes found in a Cooler file.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = RenameChromosomesConfig{};
           });

  _config = RenameChromosomesConfig{};
  auto& c = std::get<RenameChromosomesConfig>(_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a or .[ms]cool file (Cooler URI syntax supported).")
      ->required();

  sc.add_option(
      "--name-mappings",
      c.path_to_name_mappings,
      "Path to a two column TSV with pairs of chromosomes to be renamed.\n"
      "The first column should contain the original chromosome name,\n"
      "while the second column should contain the destination name to use when renaming."
  );

  sc.add_flag(
      "--add-chr-prefix",
      c.add_chr_prefix,
      "Prefix chromosome names with \"chr\".")
      ->capture_default_str();

  sc.add_flag(
      "--remove-chr-prefix",
      c.remove_chr_prefix,
      "Remove prefix \"chr\" from chromosome names.")
      ->capture_default_str();
  // clang-format on

  sc.get_option("--name-mappings")->excludes(sc.get_option("--add-chr-prefix"));
  sc.get_option("--name-mappings")->excludes(sc.get_option("--remove-chr-prefix"));
  sc.get_option("--add-chr-prefix")->excludes(sc.get_option("--remove-chr-prefix"));

  _config = std::monostate{};
}

void Cli::validate_rename_chromosomes_subcommand() const {
  assert(_cli.get_subcommand("rename-chromosomes")->parsed());

  const auto& c = std::get<RenameChromosomesConfig>(_config);

  std::vector<std::string> errors;

  if (!cooler::utils::is_cooler(c.uri) && !cooler::utils::is_multires_file(c.uri) &&
      !cooler::utils::is_scool_file(c.uri)) {
    errors.emplace_back(
        fmt::format(FMT_STRING("File \"{}\" does not appear to be a Cooler file."), c.uri));
  }

  const auto& sc = *_cli.get_subcommand("rename-chromosomes");
  if (sc.get_option("--name-mappings")->empty() && sc.get_option("--add-chr-prefix")->empty() &&
      sc.get_option("--remove-chr-prefix")->empty()) {
    errors.emplace_back(
        "please specify exactly one of --name-mappings, --add-chr-prefix, --remove-chr-prefix");
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}\n"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_rename_chromosomes_subcommand() {
  assert(_cli.get_subcommand("rename-chromosomes")->parsed());
  auto& c = std::get<RenameChromosomesConfig>(_config);

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
}

}  // namespace hictk::tools
