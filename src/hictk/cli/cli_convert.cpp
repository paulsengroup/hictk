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

void Cli::make_convert_subcommand() {
  auto& sc = *_cli.add_subcommand("convert", "Convert HiC matrices to a different format.")
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
      ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
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
      "-j,--juicer-tools-jar",
      c.juicer_tools_jar,
      "Path to juicer_tools or hic_tools JAR.")
      ->check(CLI::ExistingFile);
  sc.add_option(
      "-r,--resolutions",
      c.resolutions,
      "One or more resolution to be converted. By default all resolutions are converted.")
      ->check(CLI::PositiveNumber);
  sc.add_option(
      "--normalization-methods",
      c.normalization_methods,
      "Name of one or more normalization methods to be copied.\n"
      "By default, vectors for all known normalization methods are copied.\n")
      ->default_str("ALL");
  sc.add_flag(
      "--fail-if-norm-not-found",
      c.fail_if_normalization_method_is_not_avaliable,
      "Fail if any of the requested normalization vectors are missing.")
      ->capture_default_str();
  sc.add_option(
      "-g,--genome",
      c.genome,
      "Genome assembly name. By default this is copied from the .hic file metadata.");
  sc.add_option(
      "--juicer-tools-memory",
      c.juicer_tools_xmx,
      "Max heap size used by juicer_tools. Only used when converting from cool to hic")
      ->default_str(fmt::format(FMT_STRING("{:.0f}GB"), double(c.juicer_tools_xmx) / 1.0e9))
      ->check(CLI::PositiveNumber)
      ->transform(CLI::AsSizeValue(true));
  sc.add_option(
      "--tmpdir",
      c.tmp_dir,
      "Path where to store temporary files.");
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
      ->check(CLI::Range(std::uint32_t(2), std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_option(
      "-l,--compression-level",
      c.gzip_compression_lvl,
      "Compression level used to compress temporary files.\n"
      "Pass 0 to disable compression.")
      ->check(CLI::Range(0, 9))
      ->capture_default_str();
  sc.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files (if any).")
      ->capture_default_str();
  // clang-format on
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

    return {cooler::File(path_to_input_file.string()).bin_size()};
  }();

  std::vector<std::uint32_t> missing_resolutions;
  for (const auto res : requested_res) {
    const auto it = std::find(available_res.begin(), available_res.end(), std::int32_t(res));
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

  if ((is_cool || is_mcool) && c.juicer_tools_jar.empty()) {
    errors.emplace_back(
        fmt::format(FMT_STRING("--juicer-tools-jar is required when converting to .hic.")));
  }

  if (!c.output_format.empty()) {
    if ((is_hic && c.output_format == "hic") || (is_cool && c.output_format == "cool") ||
        (is_mcool && c.output_format == "mcool")) {
      errors.emplace_back("input and output file already are in the same format");
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

void Cli::transform_args_convert_subcommand() {
  auto& c = std::get<ConvertConfig>(_config);

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

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;

  if (c.tmp_dir.empty()) {
    c.tmp_dir = c.path_to_output.parent_path();
  }
}

}  // namespace hictk::tools
