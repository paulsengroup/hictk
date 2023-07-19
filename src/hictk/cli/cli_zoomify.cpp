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
      "Output path.");

  sc.add_flag(
      "--force",
      c.force,
      "Force overwrite existing output file(s).")
      ->capture_default_str();

  sc.add_option(
      "--resolutions",
      c.resolutions,
      "One or more resolution to be used for coarsening.")
      ->required(true);
  // clang-format on

  this->_config = std::monostate{};
}

[[nodiscard]] static std::vector<std::uint32_t> detect_duplicate_resolutions(
    const std::vector<std::uint32_t>& resolutions) {
  phmap::flat_hash_map<std::uint32_t, std::size_t> deduped_resolutions{};
  for (const auto& res : resolutions) {
    auto it = deduped_resolutions.find(res);
    if (it == deduped_resolutions.end()) {
      deduped_resolutions.emplace(res, 1);
    } else {
      it->second++;
    }
  }

  std::vector<std::uint32_t> duplicate_resolutions{};
  for (const auto& [res, count] : deduped_resolutions) {
    if (count != 1) {
      duplicate_resolutions.push_back(res);
    }
  }
  std::sort(duplicate_resolutions.begin(), duplicate_resolutions.end());
  return duplicate_resolutions;
}

static std::vector<std::uint32_t> detect_invalid_resolutions(
    const cooler::File& clr, const std::vector<std::uint32_t>& resolutions) {
  const auto base_resolution = clr.bin_size();
  std::vector<std::uint32_t> invalid_resolutions{};
  for (const auto& res : resolutions) {
    if (res % base_resolution != 0 || res < base_resolution) {
      invalid_resolutions.push_back(res);
    }
  }
  return invalid_resolutions;
}

void Cli::validate_zoomify_subcommand() const {
  assert(this->_cli.get_subcommand("zoomify")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<ZoomifyConfig>(this->_config);

  auto clr = cooler::File::open_read_only(c.input_uri);
  const auto output_path = c.output_path.empty()
                               ? std::filesystem::path(clr.path()).replace_extension(".mcool")
                               : std::filesystem::path(c.output_path);
  if (!c.force && std::filesystem::exists(output_path)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.output_path));
  }

  if (const auto dupl = detect_duplicate_resolutions(c.resolutions); !dupl.empty()) {
    errors.emplace_back(
        fmt::format(FMT_STRING("Found duplicate resolution(s):\n - {}"), fmt::join(dupl, "\n - ")));
  }

  if (const auto invalid = detect_invalid_resolutions(clr, c.resolutions); !invalid.empty()) {
    errors.emplace_back(
        fmt::format(FMT_STRING("Found the following invalid resolution(s):\n   - {}\n"
                               "Resolutions should be a multiple of the base resolution ({})."),
                    fmt::join(invalid, "\n    - "), clr.bin_size()));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n   - {}"),
                    fmt::join(errors, "\n   - ")));
  }
}

void Cli::transform_args_zoomify_subcommand() {
  auto& c = std::get<ZoomifyConfig>(this->_config);

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;

  auto clr = cooler::File::open_read_only(c.input_uri);

  if (c.output_path.empty()) {
    c.output_path = std::filesystem::path(clr.path()).replace_extension(".mcool").string();
  }

  std::sort(c.resolutions.begin(), c.resolutions.end());

  if (c.resolutions.front() != clr.bin_size()) {
    c.resolutions.insert(c.resolutions.begin(), clr.bin_size());
  }
}

}  // namespace hictk::tools
