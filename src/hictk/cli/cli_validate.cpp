//
// Created by roby on 7/13/23.
//

#include <fmt/format.h>
#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <string>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {
void Cli::make_validate_subcommand() {
  auto& sc = *this->_cli.add_subcommand("validate", "Validate .cooler and .hic files.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = ValidateConfig{};
                  });

  this->_config = ValidateConfig{};
  auto& c = std::get<ValidateConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic, .cool or .mcool file (Cooler URI syntax supported).")
      ->check(IsValidHiCFile | IsValidCoolerFile)
      ->required();

  sc.add_flag(
      "--validate-index",
      c.validate_index,
      "Validate Cooler index (may take a long time).")
      ->capture_default_str();

  sc.add_flag(
      "--quiet",
      c.validate_index,
      "Don't print anything to stdout. Success/failure is reported through exit codes")
      ->capture_default_str();
  // clang-format on

  this->_config = std::monostate{};
}
}  // namespace hictk::tools
