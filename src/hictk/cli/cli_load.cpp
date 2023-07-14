// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <string>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_load_subcommand() {
  auto& sc =
      *this->_cli
           .add_subcommand("load", "Build .cool files from interactions in various text formats.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(this->_config.index() == 0);
             this->_config = LoadConfig{};
           });

  this->_config = LoadConfig{};
  auto& c = std::get<LoadConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "chrom-sizes",
      c.path_to_chrom_sizes,
      "Path to .chrom.sizes file.")
      ->check(CLI::ExistingFile)
      ->required();

  sc.add_option(
      "bin-size",
      c.bin_size,
      "Bin size (bp).")
      ->check(CLI::PositiveNumber)
      ->required();

  sc.add_option(
      "output-uri",
      c.uri,
      "Path to output Cooler (URI syntax supported).")
      ->required();

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
      "--count-as-float",
      c.count_as_float,
      "Interactions are floats.")
      ->capture_default_str();

  sc.add_flag(
      "--assume-sorted,!--assume-unsorted",
      c.assume_sorted,
      "Assume input files are already sorted.")
      ->capture_default_str();
  sc.add_option(
      "--batch-size",
      c.batch_size,
      "Number of pixels to buffer in memory. Only used when processing unsorted interactions or pairs")
      ->capture_default_str();
  // clang-format on

  this->_config = std::monostate{};
}

void Cli::validate_load_subcommand() const {
  assert(this->_cli.get_subcommand("load")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<LoadConfig>(this->_config);

  if (!c.force && std::filesystem::exists(c.uri)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.uri));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

}  // namespace hictk::tools
