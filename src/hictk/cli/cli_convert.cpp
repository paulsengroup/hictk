// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <variant>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/multires_file.hpp"
#include "hictk/string.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/validators.hpp"

namespace hictk::tools {

void Cli::make_convert_subcommand() {
  auto& sc = *_cli.add_subcommand("convert", "Convert Hi-C files between different formats.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = ConvertConfig{};
                  });

  _config = ConvertConfig{};
  auto& c = std::get<ConvertConfig>(_config);

  // clang-format off
  sc.add_option(
      "input",
      c.path_to_input,
      "Path to the .hic, .cool or .mcool file to be converted.")
      ->check((IsValidCoolerFile | IsValidHiCFile) & (!IsValidSingleCellCoolerFile))
      ->required();
  sc.add_option(
      "output",
      c.path_to_output,
      "Output path. File extension is used to infer output format.")
      ->required();
  sc.add_option(
      "--output-fmt",
      c.output_format,
      "Output format (by default this is inferred from the output file extension).\n"
      "Should be one of:\n"
      "- cool\n"
      "- mcool\n"
      "- hic\n")
      ->check(CLI::IsMember({"cool", "mcool", "hic"}))
      ->default_str("auto");
  sc.add_option(
      "-r,--resolutions",
      c.resolutions,
      "One or more resolutions to be converted. By default all resolutions are converted.")
      ->check(CLI::PositiveNumber)
      ->transform(AsGenomicDistance);
  sc.add_option(
      "--normalization-methods",
      c.normalization_methods,
      "Name of one or more normalization methods to be copied.\n"
      "By default, vectors for all known normalization methods are copied.\n"
      "Pass NONE to avoid copying the normalization vectors.")
      ->default_str("ALL");
  sc.add_flag(
      "--fail-if-norm-not-found",
      c.fail_if_normalization_method_is_not_avaliable,
      "Fail if any of the requested normalization vectors are missing.")
      ->capture_default_str();
  sc.add_option(
      "-g,--genome",
      c.genome,
      "Genome assembly name. By default this is copied from the .hic file "
      "metadata.");
  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path where to store temporary files.")
      ->check(CLI::ExistingDirectory)
      ->capture_default_str();
  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Batch size to use when converting .[m]cool to .hic.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();
  sc.add_option(
      "-v,--verbosity",
      c.verbosity,
      "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();
  sc.add_option(
      "-t,--threads",
      c.threads,
      "Maximum number of parallel threads to spawn.\n"
      "When converting from hic to cool, only two threads will be used.")
      ->check(CLI::Range(std::uint32_t{2}, std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_option(
      "-l,--compression-lvl",
      c.compression_lvl,
      "Compression level used to compress interactions.\n"
      "Defaults to 6 and 10 for .cool and .hic files, respectively.")
      ->check(CLI::Range(std::int16_t{1}, MAX_HIC_COMPRESSION_LEVEL))
      ->capture_default_str();
  sc.add_flag(
      "--skip-all-vs-all,!--no-skip-all-vs-all",
      c.skip_all_vs_all_matrix,
      "Do not generate All vs All matrix.\n"
      "Has no effect when creating .[m]cool files.")
      ->capture_default_str();
  sc.add_option(
      "--count-type",
      c.count_type,
      "Specify the strategy used to infer count types when converting .hic "
      "files to .[m]cool format.\n"
      "Can be one of: int, float, or auto.")
      ->check(CLI::IsMember{{"auto", "int", "float"}})
      ->capture_default_str();
  sc.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files (if any).")
      ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

static void check_requested_resolutions_avail(const std::filesystem::path& path_to_input_file,
                                              const std::vector<std::uint32_t>& requested_res,
                                              std::vector<std::string>& errors) {
  const auto available_res = [&]() -> std::vector<std::uint32_t> {
    if (hic::utils::is_hic_file(path_to_input_file)) {
      return hic::utils::list_resolutions(path_to_input_file);
    }

    if (cooler::utils::is_multires_file(path_to_input_file.string())) {
      return cooler::utils::list_resolutions(path_to_input_file);
    }

    return {cooler::File(path_to_input_file.string()).resolution()};
  }();

  std::vector<std::uint32_t> missing_resolutions;
  for (const auto res : requested_res) {
    const auto it =
        std::find(available_res.begin(), available_res.end(), static_cast<std::int32_t>(res));
    if (it == available_res.end()) {
      missing_resolutions.push_back(res);
    }
  }

  if (!missing_resolutions.empty()) {
    errors.emplace_back(fmt::format(FMT_STRING("{} does not contain matrices for the following "
                                               "resolution(s): {}.\n"
                                               "Available resolutions: {}"),
                                    path_to_input_file, fmt::join(missing_resolutions, ", "),
                                    fmt::join(available_res, ", ")));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Cli::validate_convert_subcommand() const {
  const auto& c = std::get<ConvertConfig>(_config);
  std::vector<std::string> errors;

  const auto is_hic = hic::utils::is_hic_file(c.path_to_input);
  const auto is_cool = cooler::utils::is_cooler(c.path_to_input.string());
  const auto is_mcool = cooler::utils::is_multires_file(c.path_to_input.string());

  if (!is_hic && !is_cool && !is_mcool) {
    errors.emplace_back(
        fmt::format(FMT_STRING("{} is not in .hic, .cool or .mcool format"), c.path_to_input));
  }

  if (!c.output_format.empty()) {
    if ((is_hic && c.output_format == "hic") || (is_cool && c.output_format == "cool") ||
        (is_mcool && c.output_format == "mcool")) {
      errors.emplace_back("input and output file already are in the same format");
    }
  }

  const auto output_format =
      c.output_format.empty() ? infer_output_format(c.path_to_output) : c.output_format;
  if (is_cool && output_format == "hic") {
    const auto storage_mode = cooler::File(c.path_to_input.string()).attributes().storage_mode;
    if (storage_mode.has_value() && storage_mode != "symmetric-upper") {
      errors.emplace_back(fmt::format(
          FMT_STRING("converting .cool with storage-mode=\"{}\" to .hic format is not supported"),
          *storage_mode));
    }
  } else if (is_mcool && output_format == "hic") {
    const cooler::MultiResFile mclr(c.path_to_input.string());
    const auto storage_mode = mclr.open(mclr.resolutions().front()).attributes().storage_mode;
    if (storage_mode.has_value() && storage_mode != "symmetric-upper") {
      errors.emplace_back(fmt::format(
          FMT_STRING("converting .mcool with storage-mode=\"{}\" to .hic format is not supported"),
          *storage_mode));
    }
  } else if (is_hic && output_format == "cool") {
    const auto input_is_multires = MultiResFile{c.path_to_input.string()}.resolutions().size() != 1;
    if (c.resolutions.size() != 1 && input_is_multires) {
      errors.emplace_back(
          "converting multi-resolution .hic files to .cool format requires exactly one resolution "
          "to be passed through the --resolutions option");
    }
  }

  if (!c.resolutions.empty()) {
    check_requested_resolutions_avail(c.path_to_input, c.resolutions, errors);
  }

  if (!c.force && std::filesystem::exists(c.path_to_output)) {
    errors.emplace_back(fmt::format(
        FMT_STRING("Refusing to overwrite file {}. Pass --force to overwrite."), c.path_to_output));
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

// NOLINTNEXTLINE(misc-no-recursion)
[[nodiscard]] static std::string infer_assembly(const std::filesystem::path& p,
                                                std::uint32_t resolution, std::string_view format) {
  if (format == "cool") {
    auto assembly = cooler::File(p.string()).attributes().assembly;
    return !!assembly ? *assembly : "unknown";
  }
  if (format == "mcool") {
    return infer_assembly(fmt::format(FMT_STRING("{}::/resolutions/{}"), p.string(), resolution),
                          resolution, "cool");
  }
  assert(format == "hic");
  return hic::File{p.string(), resolution}.assembly();
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Cli::transform_args_convert_subcommand() {
  auto& c = std::get<ConvertConfig>(_config);
  const auto& sc = *_cli.get_subcommand("convert");

  c.input_format = infer_input_format(c.path_to_input);
  if (c.output_format.empty()) {
    c.output_format = infer_output_format(c.path_to_output);
  }

  if (c.resolutions.empty()) {
    c.resolutions = list_resolutions(c.path_to_input, c.input_format);
  }

  if (c.genome.empty()) {
    c.genome = infer_assembly(c.path_to_input, c.resolutions.back(), c.input_format);
  }

  if (c.normalization_methods.empty()) {
    if (c.input_format == "mcool") {
      c.normalization_methods =
          cooler::MultiResFile(c.path_to_input.string()).avail_normalizations("union");
    } else if (c.input_format == "cool") {
      c.normalization_methods = cooler::File(c.path_to_input.string()).avail_normalizations();
    } else {
      assert(c.input_format == "hic");
      c.normalization_methods = hic::utils::list_normalizations(c.path_to_input.string(), "union");
    }
  }

  if (sc.get_option("--tmpdir")->empty()) {
    c.tmp_dir = hictk::internal::TmpDir::default_temp_directory_path();
  }

  if (sc.get_option("--compression-lvl")->empty()) {
    c.compression_lvl =
        c.output_format == "hic" ? DEFAULT_HIC_COMPRESSION_LEVEL : DEFAULT_COOL_COMPRESSION_LEVEL;
  }

  // validate transformed args
  std::vector<std::string> errors;
  if (hictk::internal::ends_with(c.input_format, "cool") &&
      hictk::internal::ends_with(c.output_format, "cool")) {
    if (c.input_format == c.output_format) {
      errors.emplace_back(fmt::format(FMT_STRING("input is already in {} format"), c.input_format));
    } else {
      errors.emplace_back(fmt::format(FMT_STRING("converting {} -> {} is not supported"),
                                      c.input_format, c.output_format));
    }
  }

  if (hictk::internal::starts_with(c.input_format, "hic") &&
      hictk::internal::starts_with(c.output_format, "hic")) {
    errors.emplace_back("input is already in hic format");
  }

  if (hictk::internal::starts_with(c.input_format, "hic") && c.output_format == "cool" &&
      c.resolutions.size() != 1) {
    errors.emplace_back(
        "converting multi-resolution .hic files to .cool format requires specifying the resolution "
        "to be converted through the --resolutions option");
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }

  const auto try_read_from_env = sc.get_option("--verbosity")->empty();
  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
  c.verbosity = parse_hictk_verbosity_from_env(!try_read_from_env)
                    .value_or(static_cast<std::int16_t>(spdlog::level::critical) - c.verbosity);
}

}  // namespace hictk::tools
