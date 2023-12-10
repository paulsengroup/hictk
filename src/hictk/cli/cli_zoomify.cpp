// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_zoomify_subcommand() {
  auto& sc = *_cli.add_subcommand(
                      "zoomify",
                      "Convert single-resolution Cooler file to multi-resolution by coarsening.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = ZoomifyConfig{};
                  });

  _config = ZoomifyConfig{};
  auto& c = std::get<ZoomifyConfig>(_config);

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
      "One or more resolutions to be used for coarsening.");

  sc.add_flag(
      "--copy-base-resolution,!--no-copy-base-resolution",
      c.copy_base_resolution,
      "Copy the base resolution to the output file.");

  sc.add_flag(
      "--nice-steps,!--pow2-steps",
      c.nice_resolution_steps,
      "Use nice or power of two steps to automatically generate the list of resolutions.\n"
      "Example:\n"
      "Base resolution: 1000\n"
      "Pow2: 1000, 2000, 4000, 8000...\n"
      "Nice: 1000, 2000, 5000, 10000...\n")
      ->default_str("--nice-steps");

  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();

  // clang-format on

  _config = std::monostate{};
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
  assert(_cli.get_subcommand("zoomify")->parsed());

  std::vector<std::string> warnings;
  std::vector<std::string> errors;
  const auto& c = std::get<ZoomifyConfig>(_config);

  const cooler::File clr(c.input_uri);
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

  const auto* sc = _cli.get_subcommand("zoomify");
  const auto nice_or_pow2_steps_parsed =
      !sc->get_option("--nice-steps")->empty() || !sc->get_option("--pow2-steps")->empty();
  if (!c.resolutions.empty() && nice_or_pow2_steps_parsed) {
    warnings.emplace_back(
        "--nice-steps and --pow2-steps are ignored when resolutions are explicitly set with "
        "--resolutions.");
  }

  for (const auto& w : warnings) {
    SPDLOG_WARN(FMT_STRING("{}"), w);
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n   - {}"),
                    fmt::join(errors, "\n   - ")));
  }
}

static std::vector<std::uint32_t> generate_resolutions_pow2(
    std::uint32_t base_resolution, std::uint32_t upper_bound = 10'000'000) {
  assert(base_resolution != 0);
  std::vector<std::uint32_t> resolutions{base_resolution};

  for (auto res = resolutions.back(); res * 2 <= upper_bound; res = resolutions.back()) {
    resolutions.push_back(res * 2);
  }

  return resolutions;
}

static std::vector<std::uint32_t> generate_resolutions_nice(
    std::uint32_t base_resolution, std::uint32_t upper_bound = 10'000'000) {
  assert(base_resolution != 0);
  std::vector<std::uint32_t> resolutions{base_resolution};

  while (resolutions.back() * 2 <= upper_bound) {
    const auto res = resolutions.back();

    if (res * 2 > upper_bound) {
      break;
    }
    resolutions.push_back(res * 2);

    if (res * 5 > upper_bound) {
      break;
    }
    resolutions.push_back(res * 5);

    if (res * 10 > upper_bound) {
      break;
    }
    resolutions.push_back(res * 10);
  }

  return resolutions;
}

void Cli::transform_args_zoomify_subcommand() {
  auto& c = std::get<ZoomifyConfig>(_config);

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
  const cooler::File clr(c.input_uri);

  if (c.output_path.empty()) {
    c.output_path = std::filesystem::path(clr.path()).replace_extension(".mcool").string();
  }

  if (c.resolutions.empty()) {
    c.resolutions = c.nice_resolution_steps ? generate_resolutions_nice(clr.bin_size())
                                            : generate_resolutions_pow2(clr.bin_size());
  } else {
    std::sort(c.resolutions.begin(), c.resolutions.end());
  }

  if (c.resolutions.front() != clr.bin_size()) {
    c.resolutions.insert(c.resolutions.begin(), clr.bin_size());
  }
}

}  // namespace hictk::tools
