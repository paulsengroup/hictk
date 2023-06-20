// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/cli.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/common.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <regex>
#include <string>

#include "hictk/tools/config.hpp"

namespace hictk::tools {

class CoolerFileValidator : public CLI::Validator {
 public:
  CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_cooler(uri)) {
        if (hictk::cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        return "Not a valid Cooler: " + uri;
      }
      return "";
    };
  }
};

class HiCFileValidator : public CLI::Validator {
 public:
  HiCFileValidator() : Validator("HiC") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!hictk::hic::utils::is_hic_file(path)) {
        return "Not a valid .hic file: " + path;
      }
      return "";
    };
  }
};

[[nodiscard]] static std::string str_replace_all(std::string s, const std::regex& pattern,
                                                 const std::string& replacement) {
  while (std::regex_search(s, pattern)) {
    s = std::regex_replace(s, pattern, replacement);
  }
  return s;
}

class Formatter : public CLI::Formatter {
  [[nodiscard]] inline std::string make_option_opts(const CLI::Option* opt) const override {
    if (!opt->get_option_text().empty()) {
      return opt->get_option_text();
    }

    auto str_contains = [](const auto s, const auto query) {
      return s.find(query) != decltype(s)::npos;
    };

    std::string out;
    if (opt->get_type_size() != 0) {
      // Format default values so that the help string reads like: --my-option=17.0
      if (!opt->get_default_str().empty()) {
        if (internal::starts_with(opt->get_type_name(), "FLOAT")) {
          auto s = opt->get_default_str();
          if (s.find('.') == std::string::npos) {
            s += ".0";
          }
          out += fmt::format(FMT_STRING("={}"), s);
        } else {
          out += fmt::format(FMT_STRING("={}"), opt->get_default_str());
        }
      }

      // Format param domain using open/closed interval notation
      const std::regex pattern(" - ");
      if (const auto& t = opt->get_type_name(); str_contains(t, " in ")) {
        const auto p1 = t.find("[", t.find(" in "));
        const auto p2 = t.find("]", t.find(" in "));
        if (p1 != std::string::npos && p2 != std::string::npos && p2 > p1) {
          out += " " + str_replace_all(t.substr(p1, p2), pattern, ", ");
        }
      } else if (str_contains(t, "POSITIVE")) {
        out += " (0, inf)";
      } else if (str_contains(t, "NONNEGATIVE") || str_contains(t, "UINT")) {
        out += " [0, inf)";
      }

      if (opt->get_expected_max() == CLI::detail::expected_max_vector_size) {
        out += " ...";
      } else if (opt->get_expected_min() > 1) {
        out += fmt::format(FMT_STRING(" x {}"), opt->get_expected());
      }

      if (opt->get_required()) {
        out += " REQUIRED";
      }
    }
    if (!opt->get_envname().empty()) {
      out += fmt::format(FMT_STRING(" ({}: {})"), get_label("env"), opt->get_envname());
    }
    if (!opt->get_needs().empty()) {
      out += fmt::format(FMT_STRING(" {}:"), get_label("needs"));
      for (const auto* op : opt->get_needs()) {
        out += fmt::format(FMT_STRING(" {}"), op->get_name());
      }
    }
    if (!opt->get_excludes().empty()) {
      out += fmt::format(FMT_STRING(" {}:"), get_label("excludes"));
      for (const auto* op : opt->get_excludes()) {
        out += fmt::format(FMT_STRING(" {}"), op->get_name());
      }
    }

    return out;
  }
};

inline const auto IsValidCoolerFile = CoolerFileValidator();
inline const auto IsValidHiCFile = HiCFileValidator();

// clang-format off
inline const auto ParseHiCMatrixType = CLI::CheckedTransformer(
    std::map<std::string, hictk::hic::MatrixType>{
        {"observed", hictk::hic::MatrixType::observed},
        {"oe", hictk::hic::MatrixType::oe},
        {"expected", hictk::hic::MatrixType::expected}},
    CLI::ignore_case);

inline const auto ParseHiCNormalization = CLI::CheckedTransformer(
    std::map<std::string, hictk::hic::NormalizationMethod>{
        {"NONE", hictk::hic::NormalizationMethod::NONE},
        {"VC", hictk::hic::NormalizationMethod::VC},
        {"VC_SQRT", hictk::hic::NormalizationMethod::VC_SQRT},
        {"KR", hictk::hic::NormalizationMethod::KR},
        {"SCALE", hictk::hic::NormalizationMethod::SCALE},
        {"INTER_VC", hictk::hic::NormalizationMethod::INTER_VC},
        {"INTER_KR", hictk::hic::NormalizationMethod::INTER_KR},
        {"INTER_SCALE", hictk::hic::NormalizationMethod::INTER_SCALE},
        {"GW_VC", hictk::hic::NormalizationMethod::GW_VC},
        {"GW_KR", hictk::hic::NormalizationMethod::GW_KR},
        {"GW_SCALE", hictk::hic::NormalizationMethod::GW_SCALE}},
    CLI::ignore_case);

inline const auto ParseHiCMatrixUnit = CLI::CheckedTransformer(
    std::map<std::string, hictk::hic::MatrixUnit>{
        {"BP", hictk::hic::MatrixUnit::BP},
        {"FRAG", hictk::hic::MatrixUnit::FRAG}},
    CLI::ignore_case);
// clang-format on

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return this->_subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(this->get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  try {
    this->_cli.name(this->_exec_name);
    this->_cli.parse(this->_argc, this->_argv);

    if (this->_cli.get_subcommand("convert")->parsed()) {
      this->_subcommand = subcommand::convert;
    } else if (this->_cli.get_subcommand("dump")->parsed()) {
      this->_subcommand = subcommand::dump;
    } else if (this->_cli.get_subcommand("load")->parsed()) {
      this->_subcommand = subcommand::load;
    } else if (this->_cli.get_subcommand("merge")->parsed()) {
      this->_subcommand = subcommand::merge;
    } else {
      this->_subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    this->_exit_code = this->_cli.exit(e);
    return this->_config;
  } catch (const std::exception& e) {
    this->_exit_code = 1;
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
            "message, please file an issue on GitHub"),
        e.what()));

  } catch (...) {
    this->_exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  this->validate();
  this->transform_args();

  this->_exit_code = 0;
  return this->_config;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  switch (s) {
    case convert:
      return "convert";
    case dump:
      return "dump";
    case load:
      return "load";
    case merge:
      return "merge";
    default:
      assert(s == help);
      return "--help";
  }
}

void Cli::make_convert_subcommand() {
  auto& sc = *this->_cli.add_subcommand("convert", "Convert HiC matrices to a different format.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = ConvertConfig{};
                  });

  this->_config = ConvertConfig{};
  auto& c = std::get<ConvertConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "hic",
      c.input_hic, "Path to the .hic file to be converted.")
      ->check(CLI::ExistingFile)
      ->required();
  sc.add_option(
      "-o,--output",
      c.output_cooler,
      "Path where to store the output cooler.");
  sc.add_option(
      "-r,--resolutions",
      c.resolutions,
      "One or more resolution to be converted. By default all resolutions are converted.")
      ->check(CLI::PositiveNumber);
  sc.add_option(
       "--normalization-methods",
       c.normalization_methods_str,
       fmt::format(FMT_STRING("Name of one or more normalization methods to be copied.\n"
                              "By default vectors for all known normalization methods are copied.\n"
                              "Supported methods:\n - {}"),
                   fmt::join(hic::NORMALIZATION_METHODS, "\n - ")))
      ->default_str("ALL");
  sc.add_flag(
      "--fail-if-norm-not-found",
      c.fail_if_normalization_method_is_not_avaliable,
      "Fail if any of the requested normalization vectors are missing.")
      ->capture_default_str();
  sc.add_option("-g,--genome", c.genome,
               "Genome assembly name. By default this is copied from the .hic file metadata.");
  sc.add_flag("-q,--quiet", c.quiet, "Suppress console output.")->capture_default_str();
  sc.add_option("-v,--verbosity", c.verbosity, "Set verbosity of output to the console.")
      ->check(CLI::Range(1, 4))
      ->excludes("--quiet")
      ->capture_default_str();
  sc.add_flag("-f,--force", c.force, "Overwrite existing files (if any).")->capture_default_str();
  // clang-format on
}

void Cli::make_dump_subcommand() {
  auto& sc = *this->_cli.add_subcommand("dump", "Dump Cooler data to stdout.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = DumpConfig{};
                  });

  this->_config = DumpConfig{};
  auto& c = std::get<DumpConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic, .cool or .mcool file (Cooler URI syntax supported).")
      ->check(IsValidHiCFile | IsValidCoolerFile)
      ->required();

  sc.add_option(
      "--resolution",
      c.resolution,
      "HiC matrix resolution (ignored when file is not in .hic format).")
      ->check(CLI::NonNegativeNumber);

  sc.add_option(
      "--matrix-type",
      c.matrix_type,
      "Matrix type (ignored when file is not in .hic format).")
      ->transform(ParseHiCMatrixType)
      ->default_str("observed");

  sc.add_option(
      "--matrix-unit",
      c.matrix_unit,
      "Matrix unit (ignored when file is not in .hic format).")
      ->transform(ParseHiCMatrixUnit)
      ->default_str("BP");

  sc.add_option(
      "-t,--table",
      c.table,
      "Name of the table to dump.\n")
      ->check(CLI::IsMember({"chroms", "bins", "pixels"}))
      ->capture_default_str();

  sc.add_option(
      "-r,--range",
      c.range1,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--range2",
      c.range2,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--query-file",
      c.query_file,
      "Path to a BEDPE file with the list of coordinates to be fetched (pass - to read queries from stdin).")
      ->check(CLI::ExistingFile | CLI::IsMember({"-"}))
      ->capture_default_str();

  sc.add_option(
      "-b,--balance",
      c.normalization,
      "Balance interactions using the given method.")
      ->capture_default_str();

  sc.add_flag(
      "--join,!--no-join",
      c.join,
      "Output pixels in BG2 format.")
      ->capture_default_str();

  sc.add_option(
      "--weight-type",
      c.weight_type,
      "How balancing weights should be applied to raw interactions (ignored when file is in .hic format).")
      ->check(CLI::IsMember({"infer", "divisive", "multiplicative"}))
      ->capture_default_str();

  // clang-format on

  sc.get_option("--query-file")->excludes(sc.get_option("--range"));
  sc.get_option("--query-file")->excludes(sc.get_option("--range2"));

  this->_config = std::monostate{};
}

void Cli::make_load_subcommand() {
  auto& sc = *this->_cli.add_subcommand("load", "Build .cool files from interactions in BG2/COO.")
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
      ->check(CLI::IsMember({"bg2", "coo"}))
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
      "--assume-assume_sorted,!--no-assume-sorted",
      c.assume_sorted,
      "Assume input files are already assume_sorted.")
      ->capture_default_str();
  // clang-format on

  this->_config = std::monostate{};
}

void Cli::make_merge_subcommand() {
  auto& sc = *this->_cli.add_subcommand("merge", "Merge coolers.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = MergeConfig{};
                  });

  this->_config = MergeConfig{};
  auto& c = std::get<MergeConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "input-coolers",
      c.input_uris,
      "Path to two or more Cooler files to be merged (URI syntax supported).")
      ->check(IsValidCoolerFile)
      ->expected(2, std::numeric_limits<int>::max())
      ->required();

  sc.add_option(
      "-o,--output-cooler",
      c.output_uri,
      "Output Cooler (URI syntax supported).\n"
      "When not specified, merged interactions will be printed to stdout.");

  sc.add_flag(
      "-f,--force",
      c.force,
      "Force overwrite output cooler.")
      ->capture_default_str();

  sc.add_option(
      "--chunk-size",
      c.chunk_size,
      "Number of pixels to store in memory before writing to disk.")
      ->capture_default_str();

  // clang-format on

  this->_config = std::monostate{};
}

void Cli::make_cli() {
  this->_cli.name(this->_exec_name);
  this->_cli.description("Coolerpp tools.");
  this->_cli.set_version_flag("-V,--version", "hictk::cooler-tools-0.0.1");
  this->_cli.require_subcommand(1);

  this->make_convert_subcommand();
  this->make_dump_subcommand();
  this->make_load_subcommand();
  this->make_merge_subcommand();
}

static void check_requested_resolutions_avail(const std::filesystem::path& path_to_hic,
                                              const std::vector<std::uint32_t>& requested_res,
                                              std::vector<std::string>& errors) {
  const auto available_res = hic::utils::list_resolutions(path_to_hic);
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
                                    path_to_hic, fmt::join(missing_resolutions, ", "),
                                    fmt::join(available_res, ", ")));
  }
}

static void check_normalization_methods(const std::vector<std::string>& norm_methods,
                                        std::vector<std::string>& errors) {
  assert(!norm_methods.empty());
  if (norm_methods.size() > 1 && norm_methods.front() != "ALL") {
    for (const auto& n : norm_methods) {
      try {
        std::ignore = hic::ParseNormStr(n);
      } catch (...) {
        errors.emplace_back(
            fmt::format(FMT_STRING("\"{}\" is not a known normalization method"), n));
      }
    }
  }
}

void Cli::validate_convert_subcommand() const {
  auto& c = std::get<ConvertConfig>(this->_config);
  std::vector<std::string> errors;

  if (!hic::utils::is_hic_file(c.input_hic)) {
    errors.emplace_back(
        fmt::format(FMT_STRING("{} does not appear to be a valid .hic file"), c.input_hic));
  }

  if (!c.resolutions.empty()) {
    check_requested_resolutions_avail(c.input_hic, c.resolutions, errors);
  }

  check_normalization_methods(c.normalization_methods_str, errors);

  if (!errors.empty()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error(s) where encountered while validating CLI arguments:\n - {}"),
        fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_dump_subcommand() const {
  assert(this->_cli.get_subcommand("dump")->parsed());

  [[maybe_unused]] std::vector<std::string> warnings;  // TODO issue warnings
  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }

  const auto is_hic = hic::utils::is_hic_file(c.uri);
  const auto is_cooler = cooler::utils::is_cooler(c.uri);
  const auto is_mcooler = cooler::utils::is_multires_file(c.uri);

  if (is_hic && c.resolution == 0) {
    errors.push_back("--resolution is mandatory when file is in .hic format.");
  }

  if ((is_cooler || is_mcooler) && c.resolution != 0) {
    warnings.push_back("--resolution is ignored when file is in .cool or .mcool format.");
  }

  if (is_hic && c.weight_type == "infer") {
    warnings.push_back("--weight-type is ignored when file is in .hic format.");
  }

  const auto matrix_type_parsed =
      this->_cli.get_subcommand("dump")->get_option("--matrix-type")->empty();
  const auto matrix_unit_parsed =
      this->_cli.get_subcommand("dump")->get_option("--matrix-unit")->empty();

  if (!is_hic && (matrix_type_parsed || matrix_unit_parsed)) {
    warnings.push_back(
        "--matrix-type and --matrix-unit are ignored when file is not in .hic format.");
  }

  if (is_hic && c.matrix_unit == hic::MatrixUnit::FRAG) {
    errors.push_back("--matrix-type=FRAG is not yet supported.");
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}\n"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_load_subcommand() const {
  assert(this->_cli.get_subcommand("load")->parsed());
  /*
  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
  */
}

void Cli::validate_merge_subcommand() const {
  assert(this->_cli.get_subcommand("merge")->parsed());
  /*
  std::vector<std::string> errors;
  const auto& c = std::get<MergeConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
  */
}

void Cli::validate() const {
  switch (this->_subcommand) {
    case convert:
      this->validate_convert_subcommand();
      break;
    case dump:
      this->validate_dump_subcommand();
      break;
    case load:
      this->validate_load_subcommand();
      break;
    case merge:
      this->validate_merge_subcommand();
      break;
    case help:
      break;
  }
}

[[nodiscard]] static std::vector<hic::NormalizationMethod> generate_norm_vect(
    const std::vector<std::string>& norms_str) {
  assert(!norms_str.empty());
  if (norms_str.size() == 1 && norms_str.front() == "ALL") {
    return {hic::NORMALIZATION_METHODS.begin(), hic::NORMALIZATION_METHODS.end()};
  }

  std::vector<hic::NormalizationMethod> norms(norms_str.size());
  std::transform(norms_str.begin(), norms_str.end(), norms.begin(),
                 [&](const auto& s) { return hic::ParseNormStr(s); });
  return norms;
}

void Cli::transform_args_convert_subcommand() {
  auto& c = std::get<ConvertConfig>(this->_config);
  if (c.resolutions.empty()) {
    c.resolutions = hic::utils::list_resolutions(c.input_hic);
  }

  c.normalization_methods = generate_norm_vect(c.normalization_methods_str);

  if (c.genome.empty()) {
    const hic::HiCFile f(c.input_hic, c.resolutions.front());
    c.genome = f.assembly();
  }

  if (c.output_cooler.empty()) {
    const auto* extension = c.resolutions.size() > 1 ? ".mcool" : ".cool";
    c.output_cooler = c.input_hic;
    c.output_cooler.replace_extension(extension);
  }

  if (c.norm_dset_names.empty()) {
    c.norm_dset_names.resize(c.normalization_methods.size());
    std::transform(c.normalization_methods.begin(), c.normalization_methods.end(),
                   c.norm_dset_names.begin(), [](const auto norm) { return fmt::to_string(norm); });
  }

  if (c.block_cache_size != 0) {
    c.block_cache_size =
        static_cast<std::size_t>(double(c.block_cache_size) * 0.85) / sizeof(hic::SerializedPixel);
  }

  if (c.quiet) {
    c.verbosity = spdlog::level::err;
  } else {
    // in spdlog, high numbers correspond to low log levels
    assert(c.verbosity > 0 && c.verbosity < 5);
    c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;
  }
}

void Cli::transform_args() {
  switch (this->_subcommand) {
    case convert:
      this->transform_args_convert_subcommand();
      break;
    case dump:
      // this->transform_args_dump_subcommand();
      break;
    case load:
      // this->transform_args_load_subcommand();
      break;
    case merge:
      // this->transform_args_merge_subcommand();
      break;
    case help:
      break;
  }
}

}  // namespace hictk::tools
