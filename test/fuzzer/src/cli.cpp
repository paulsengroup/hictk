// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/fuzzer/cli.hpp"

#ifdef _WIN32
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>

#include <climits>

#else
#include <unistd.h>
#endif

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <array>
#include <cassert>
#include <exception>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

#include "hictk/balancing/methods.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/fuzzer/config.hpp"
#include "hictk/multires_file.hpp"
#include "hictk/version.hpp"

namespace hictk::fuzzer {

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  try {
    _cli.name(_exec_name);
    _cli.parse(_argc, _argv);

    if (_cli.get_subcommand("fuzz")->parsed()) {
      _subcommand = subcommand::fuzz;
    } else if (_cli.get_subcommand("launch-worker")->parsed()) {
      _subcommand = subcommand::launch_worker;
    } else {
      _subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    _exit_code = _cli.exit(e);
    return _config;
  } catch (const std::exception& e) {
    _exit_code = 1;
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
            "message, please file an issue on GitHub"),
        e.what()));

  } catch (...) {
    _exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  validate_args();
  transform_args();

  _exit_code = 0;
  return _config;
}

int Cli::exit(const CLI::ParseError& e) const { return _cli.exit(e); }
int Cli::exit() const noexcept { return _exit_code; }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  using sc = subcommand;
  switch (s) {
    case sc::fuzz:
      return "fuzz";
    case sc::launch_worker:
      return "launch-worker";
    default:
      assert(s == sc::help);
      return "--help";
  }
}

static void add_common_args(CLI::App& sc, Config& c) {
  sc.add_option("test-uri", c.test_uri,
                "Path to the .hic, .cool or .mcool file to be used as test file.")
      ->check(IsValidHiCFile | IsValidCoolerFile | IsValidMultiresCoolerFile)
      ->required();
  sc.add_option("reference-uri", c.reference_uri,
                "Path to the .cool or .mcool file to be used as reference file.")
      ->check(IsValidCoolerFile | IsValidMultiresCoolerFile)
      ->required();
  sc.add_option("--resolution", c.resolution,
                "Matrix resolution.\n"
                "Required when either test-uri or reference-uri are multi-resolution files.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_option("--1d-to-2d-query-ratio", c._1d_to_2d_query_ratio,
                "Ratio of 1D to 2D queries. Use 0 or 1 to only test 1D or 2D queries.")
      ->check(CLI::Bound(0.0, 1.0))
      ->capture_default_str();
  sc.add_option("--duration", c.duration, "Test duration in seconds.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();
  sc.add_option("--format", c.query_format, "Format used to fetch pixels.")
      ->check(CLI::IsMember{{"dense", "df", "sparse"}})
      ->capture_default_str();
  sc.add_option("--query-length-avg", c.query_length_avg, "Average query size.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_option("--query-length-std", c.query_length_std, "Query size standard deviation.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  sc.add_option("--normalization", c.normalization, "Name of the dataset to use for balancing.")
      ->capture_default_str();
  sc.add_flag("--join", c.join,
              "Fetch pixels in BG2 format.\n"
              "Ignored when --format is not df.")
      ->capture_default_str();
  sc.add_option("--diagonal-band-width", c.diagonal_band_width, "Diagonal band width.")
      ->capture_default_str();
  sc.add_option("--seed", c.seed, "Seed used for PRNG.")->capture_default_str();
  sc.add_option("-v,--verbosity", c.verbosity, "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->capture_default_str();
}

void Cli::make_fuzz_subcommand() {
  auto& sc = *_cli.add_subcommand("fuzz", "Run the fuzz test suite.")->fallthrough();

  auto& c = _config;
  add_common_args(sc, _config);
  sc.add_option("--nproc", c.nproc, "Number of test processes to run in parallel.")
      ->check(CLI::Bound(1U, std::thread::hardware_concurrency()))
      ->capture_default_str();
  sc.add_flag("--suppress-py-warnings,!--show-py-warnings", c.suppress_python_warnings,
              "Hide/show python warnings.")
      ->capture_default_str();
}

void Cli::make_launch_worker_subcommand() {
  auto& sc = *_cli.add_subcommand("launch-worker", "Lunch one instance of the fuzzer process.")
                  ->group("")
                  ->fallthrough();

  sc.add_option("--task-id", _config.task_id, "Task ID.")->check(CLI::PositiveNumber)->required();

  add_common_args(sc, _config);
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description("Fuzzer for hictk.");
  _cli.set_version_flag("-V,--version", "0.0.1");
  _cli.require_subcommand(1);

  make_fuzz_subcommand();
  make_launch_worker_subcommand();
}

static void validate_resolution(const std::filesystem::path& uri, std::uint32_t expected_resolution,
                                std::vector<std::string>& errors) {
  if (expected_resolution != 0) {
    if (hictk::cooler::utils::is_cooler(uri.string())) {
      const auto found_resolution = cooler::File(uri.string()).resolution();
      if (found_resolution != expected_resolution) {
        errors.emplace_back(fmt::format(FMT_STRING("expected resolution {}, found {} at URI {}"),
                                        expected_resolution, found_resolution, uri));
      }
      return;
    }

    const auto resolutions = MultiResFile(uri).resolutions();
    if (std::find(resolutions.begin(), resolutions.end(), expected_resolution) ==
        resolutions.end()) {
      errors.emplace_back(
          fmt::format(FMT_STRING("file at URI {} does not contain interactions for resolution {}"),
                      uri, expected_resolution));
    }
    return;
  }

  if (!hictk::cooler::utils::is_cooler(uri.string())) {
    errors.emplace_back(
        fmt::format(FMT_STRING("URI {} does not point to a single-resolution Cooler file: "
                               "--resolution is required when providing .hic or .mcool files"),
                    uri));
  }
}

static void validate_normalization(const std::filesystem::path& uri, std::uint32_t resolution,
                                   std::string_view normalization,
                                   std::vector<std::string>& errors) {
  if (normalization == "NONE") {
    return;
  }

  const auto avail_normalizations = File(uri.string(), resolution).avail_normalizations();
  if (std::find(avail_normalizations.begin(), avail_normalizations.end(), normalization) ==
      avail_normalizations.end()) {
    errors.emplace_back(fmt::format(
        FMT_STRING("file {} does not contain \"{}\" balancing weights at resolution {}"), uri,
        normalization, resolution));
  }
}

static void validate_common_args(const Config& c) {
  std::vector<std::string> errors;

  validate_resolution(c.test_uri, c.resolution, errors);
  validate_resolution(c.reference_uri, c.resolution, errors);

  if (errors.empty()) {
    validate_normalization(c.test_uri, c.resolution, c.normalization, errors);
    validate_normalization(c.reference_uri, c.resolution, c.normalization, errors);
  }

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_fuzz_subcommand() const {
  const auto& c = _config;

  validate_common_args(c);

  std::vector<std::string> warnings;

  if (c.diagonal_band_width.has_value() && c.query_format != "df") {
    warnings.emplace_back("--diagonal-band-width is ignored when --format is not \"df\"");
  }

  if (c.join && c.query_format != "df") {
    warnings.emplace_back("--join is ignored when --format is not \"df\"");
  }

  for (const auto& w : warnings) {
    SPDLOG_WARN(w);
  }
}

void Cli::validate_launch_worker_subcommand() const { validate_common_args(_config); }

void Cli::validate_args() const {
  using sc = subcommand;
  switch (_subcommand) {
    case sc::fuzz:
      validate_fuzz_subcommand();
      break;
    case sc::launch_worker:
      validate_launch_worker_subcommand();
      break;
    case sc::help:
      break;
  }
}

static std::string get_path_to_executable() {
#ifdef _WIN32
  std::string path(MAX_PATH, '\0');
  if (GetModuleFileNameA(NULL, path.data(), path.size())) {
    return std::string{path.c_str()};
  }

#elif defined(__APPLE__)
  std::string path(PATH_MAX, '\0');
  std::uint32_t count = PATH_MAX;
  if (!_NSGetExecutablePath(path.data(), &count)) {
    return path.substr(0, count);
  }

#else
  std::string path(PATH_MAX, '\0');
  const auto count = readlink("/proc/self/exe", path.data(), path.size());
  if (count != -1) {
    return path.substr(0, static_cast<std::size_t>(count));
  }
#endif
  throw std::runtime_error("unable to generate the path to hictk_fuzzer.");
}

void Cli::transform_args_fuzz_subcommand() {
  _config.exec = get_path_to_executable();

  if (!_config.seed.has_value()) {
    _config.seed = 0;
    // NOLINTNEXTLINE(*-reinterpret-cast)
    auto* ptr = reinterpret_cast<std::uint32_t*>(&_config.seed.value());
    *ptr++ = std::random_device{}();  // NOLINT(*-pointer-arithmetic)
    *ptr = std::random_device{}();
  }
}

void Cli::transform_args_launch_worker_subcommand() {
  // in spdlog, high numbers correspond to low log levels
  assert(_config.verbosity > 0 &&
         _config.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
  _config.verbosity = static_cast<std::int16_t>(spdlog::level::critical) - _config.verbosity;
}

void Cli::transform_args() {
  using sc = subcommand;
  switch (_subcommand) {
    case sc::fuzz:
      transform_args_fuzz_subcommand();
      break;
    case sc::launch_worker:
      transform_args_launch_worker_subcommand();
      break;
    case sc::help:
      break;
  }
}

}  // namespace hictk::fuzzer
