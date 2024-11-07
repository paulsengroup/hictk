// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
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
#include "hictk/tmpdir.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_zoomify_subcommand() {
  auto& sc =
      *_cli.add_subcommand(
               "zoomify",
               "Convert single-resolution Cooler and .hic files to multi-resolution by coarsening.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
             assert(_config.index() == 0);
             _config = ZoomifyConfig{};
           });

  _config = ZoomifyConfig{};
  auto& c = std::get<ZoomifyConfig>(_config);

  // clang-format off
  sc.add_option(
      "cooler/hic",
       c.path_to_input,
      "Path to a .cool or .hic file (Cooler URI syntax supported).")
      ->check(IsValidCoolerFile | IsValidHiCFile)
      ->required();

  sc.add_option(
      "[m]cool/hic",
      c.path_to_output,
      "Output path.\n"
      "When zoomifying Cooler files, providing a single resolution through\n"
      "--resolutions and specifying --no-copy-base-resolution, the output file\n"
      "will be in .cool format.")
      ->required();

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
      "-l,--compression-lvl",
      c.compression_lvl,
      "Compression level used to compress interactions.\n"
      "Defaults to 6 and 10 for .mcool and .hic files, respectively.")
      ->check(CLI::Bound(std::int16_t{1}, MAX_HIC_COMPRESSION_LEVEL))
      ->capture_default_str();

  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.\n"
      "When zoomifying interactions from a .cool file, only a single thread will be used.")
            ->check(CLI::Range(std::uint32_t{1}, std::thread::hardware_concurrency()))
            ->capture_default_str();

  sc.add_option(
      "--chunk-size",
      c.batch_size,
      "Number of pixels to buffer in memory.\n"
      "Only used when zoomifying .hic files.")
      ->capture_default_str();

  sc.add_flag(
      "--skip-all-vs-all,!--no-skip-all-vs-all",
      c.skip_all_vs_all_matrix,
      "Do not generate All vs All matrix.\n"
      "Has no effect when zoomifying .cool files.")
      ->capture_default_str();

  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path to a folder where to store temporary data.")
      ->check(CLI::ExistingDirectory)
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
    std::uint32_t base_resolution, const std::vector<std::uint32_t>& resolutions) {
  std::vector<std::uint32_t> invalid_resolutions{};
  for (const auto& res : resolutions) {
    if (res % base_resolution != 0 || res < base_resolution) {
      invalid_resolutions.push_back(res);
    }
  }
  return invalid_resolutions;
}

[[nodiscard]] static std::uint32_t detect_base_resolution(std::string_view path,
                                                          std::string_view format) {
  if (format == "cool") {
    return cooler::File(path).resolution();
  }

  assert(format == "hic");
  return hic::utils::list_resolutions(std::string{path}, true).front();
}

void Cli::validate_zoomify_subcommand() const {
  assert(_cli.get_subcommand("zoomify")->parsed());

  std::vector<std::string> warnings;
  std::vector<std::string> errors;
  const auto& c = std::get<ZoomifyConfig>(_config);

  if (!c.force && std::filesystem::exists(c.path_to_output)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.path_to_output));
  }

  const auto input_format = infer_input_format(c.path_to_input);
  const auto output_format = infer_output_format(c.path_to_output);
  if ((input_format == "hic" && output_format != "hic") ||
      (input_format != "hic" && output_format == "hic")) {
    errors.emplace_back(
        fmt::format(FMT_STRING("Zoomifying a .{} file to produce .{} file is not supported."),
                    input_format, output_format));
  }

  if (input_format == "cool") {
    if (const auto storage_mode = cooler::File(c.path_to_input.string()).attributes().storage_mode;
        storage_mode.has_value() && storage_mode != "symmetric-upper") {
      errors.emplace_back(fmt::format(
          FMT_STRING("Zoomifying .cool files with storage-mode=\"{}\" is not supported."),
          *storage_mode));
    }
  }

  const auto base_resolution = detect_base_resolution(c.path_to_input.string(), input_format);

  if (base_resolution == 0) {  // Variable bin size
    errors.clear();
    warnings.clear();
    errors.emplace_back("Zoomifying files with variable bin size is currently not supported.");
  } else {
    if (const auto dupl = detect_duplicate_resolutions(c.resolutions); !dupl.empty()) {
      errors.emplace_back(fmt::format(FMT_STRING("Found duplicate resolution(s):\n - {}"),
                                      fmt::join(dupl, "\n - ")));
    }

    if (const auto invalid = detect_invalid_resolutions(base_resolution, c.resolutions);
        !invalid.empty()) {
      errors.emplace_back(
          fmt::format(FMT_STRING("Found the following invalid resolution(s):\n   - {}\n"
                                 "Resolutions should be a multiple of the base resolution ({})."),
                      fmt::join(invalid, "\n    - "), base_resolution));
    }

    const auto* sc = _cli.get_subcommand("zoomify");
    const auto nice_or_pow2_steps_parsed =
        !sc->get_option("--nice-steps")->empty() || !sc->get_option("--pow2-steps")->empty();
    if (!c.resolutions.empty() && nice_or_pow2_steps_parsed) {
      warnings.emplace_back(
          "--nice-steps and --pow2-steps are ignored when resolutions are explicitly set with "
          "--resolutions.");
    }
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

static std::vector<std::uint32_t>
generate_resolutions_pow2(  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::uint32_t base_resolution, std::uint32_t upper_bound = 10'000'000) {
  assert(base_resolution != 0);
  std::vector<std::uint32_t> resolutions{base_resolution};

  for (auto res = resolutions.back(); res * 2 <= upper_bound; res = resolutions.back()) {
    resolutions.push_back(res * 2);
  }

  return resolutions;
}

static std::vector<std::uint32_t>
generate_resolutions_nice(  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::uint32_t base_resolution, std::uint32_t upper_bound = 10'000'000) {
  assert(base_resolution != 0);
  std::vector<std::uint32_t> resolutions{base_resolution};

  // NOLINTBEGIN(*-avoid-magic-numbers)
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
  // NOLINTEND(*-avoid-magic-numbers)

  return resolutions;
}

void Cli::transform_args_zoomify_subcommand() {
  auto& c = std::get<ZoomifyConfig>(_config);
  const auto& sc = *_cli.get_subcommand("zoomify");

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
  c.verbosity = static_cast<std::int16_t>(spdlog::level::critical) - c.verbosity;

  c.input_format = infer_input_format(c.path_to_input);
  c.output_format = infer_output_format(c.path_to_output);

  const auto base_resolution = detect_base_resolution(c.path_to_input.string(), c.input_format);

  if (c.resolutions.empty()) {
    c.resolutions = c.nice_resolution_steps ? generate_resolutions_nice(base_resolution)
                                            : generate_resolutions_pow2(base_resolution);
  } else {
    std::sort(c.resolutions.begin(), c.resolutions.end());
  }

  if (c.output_format == "cool" && c.resolutions.front() != base_resolution) {
    c.resolutions.insert(c.resolutions.begin(), base_resolution);
  }

  if (sc.get_option("--compression-lvl")->empty()) {
    c.compression_lvl =
        c.output_format == "hic" ? DEFAULT_HIC_COMPRESSION_LEVEL : DEFAULT_COOL_COMPRESSION_LEVEL;
  }

  if (sc.get_option("--tmpdir")->empty()) {
    c.tmp_dir = hictk::internal::TmpDir::default_temp_directory_path();
  }
}

}  // namespace hictk::tools
