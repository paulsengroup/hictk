// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <CLI/CLI.hpp>
#include <cassert>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_zoomify_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("zoomify",
                                  "Convert single-resolution cooler file to multi-resolution "
                                  "cooler file by coarsening.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = ZoomifyConfig{};
                  });

  this->_config = ZoomifyConfig{};
  auto& c = std::get<ZoomifyConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "cooler",
       c.input_uri,
      "Path to a .cool file (Cooler URI syntax supported).")
      ->check(IsValidCoolerFile)
      ->required();

  sc.add_option(
      "mcool",
      c.output_path,
      "Output path.")
      ->required();

  sc.add_option(
      "--resolutions",
      c.resolutions,
      "One or more resolution to be used for coarsening.")
      ->required(true);
  // clang-format on

  this->_config = std::monostate{};
}

}  // namespace hictk::tools
