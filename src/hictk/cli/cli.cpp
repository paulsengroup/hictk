// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/cli.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/common.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstdint>
#include <regex>
#include <string>

#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
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

class MultiresCoolerFileValidator : public CLI::Validator {
 public:
  MultiresCoolerFileValidator() : Validator("Multires-cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_multires_file(uri)) {
        return "Not a valid multi-resolution cooler: " + uri;
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
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
        const auto p1 = t.find('[', t.find(" in "));
        const auto p2 = t.find(']', t.find(" in "));
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

// clang-format off
inline const auto IsValidCoolerFile = CoolerFileValidator();                  // NOLINT(cert-err58-cpp)
inline const auto IsValidMultiresCoolerFile = MultiresCoolerFileValidator();  // NOLINT(cert-err58-cpp)
inline const auto IsValidHiCFile = HiCFileValidator();                        // NOLINT(cert-err58-cpp)
// clang-format on

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
inline const auto ParseHiCMatrixType = CLI::CheckedTransformer(
    std::map<std::string, hictk::hic::MatrixType>{
        {"observed", hictk::hic::MatrixType::observed},
        {"oe", hictk::hic::MatrixType::oe},
        {"expected", hictk::hic::MatrixType::expected}},
    CLI::ignore_case);

// NOLINTNEXTLINE(cert-err58-cpp)
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

// NOLINTNEXTLINE(cert-err58-cpp)
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
    } else if (this->_cli.get_subcommand("validate")->parsed()) {
      this->_subcommand = subcommand::validate;
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
  this->validate_args();
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
    case validate:
      return "validate";
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
       c.normalization_methods_str,
       fmt::format(FMT_STRING("Name of one or more normalization methods to be copied.\n"
                              "By default vectors for all known normalization methods are copied.\n"
                              "Supported methods:\n"
                              "- weights\n"
                              "- {}"),
                   fmt::join(hic::NORMALIZATION_METHODS, "\n - ")))
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
      "--read-cache-size",
      c.block_cache_size,
      "Maximum size of the in-memory read cache. Not used when converting to .hic")
      ->default_str("auto")
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
      "-p,--processes",
      c.processes,
      "Maximum number of parallel processes to spawn.\n"
      "When converting from hic to cool, only two processes will be used.")
      ->check(CLI::Range(2, 1024))
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

void Cli::make_cli() {
  this->_cli.name(this->_exec_name);
  this->_cli.description("Coolerpp tools.");
  this->_cli.set_version_flag("-V,--version", "hictk::cooler-tools-0.0.1");
  this->_cli.require_subcommand(1);

  this->make_convert_subcommand();
  this->make_dump_subcommand();
  this->make_load_subcommand();
  this->make_merge_subcommand();
  this->make_validate_subcommand();
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

    return {cooler::File::open_read_only(path_to_input_file.string()).bin_size()};
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
  const auto& c = std::get<ConvertConfig>(this->_config);
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

  [[maybe_unused]] std::vector<std::string> warnings;
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

  if (is_hic && c.resolution == 0 && c.table != "chroms") {
    errors.emplace_back("--resolution is mandatory when file is in .hic format.");
  }

  const auto resolution_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--resolution")->empty();

  if ((is_cooler || is_mcooler) && resolution_parsed) {
    warnings.emplace_back("--resolution is ignored when file is in .cool or .mcool format.");
  }

  const auto weight_type_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--weight-type")->empty();

  if (is_hic && weight_type_parsed) {
    warnings.emplace_back("--weight-type is ignored when file is in .hic format.");
  }

  const auto matrix_type_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--matrix-type")->empty();
  const auto matrix_unit_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--matrix-unit")->empty();

  if (!is_hic && (matrix_type_parsed || matrix_unit_parsed)) {
    warnings.emplace_back(
        "--matrix-type and --matrix-unit are ignored when input file is not in .hic format.");
  }

  if (is_hic && c.matrix_unit == hic::MatrixUnit::FRAG) {
    errors.emplace_back("--matrix-type=FRAG is not yet supported.");
  }

  for (const auto& w : warnings) {
    spdlog::warn(FMT_STRING("{}"), w);
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

void Cli::validate_args() const {
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
    case validate:
      [[fallthrough]];
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

[[nodiscard]] static std::string infer_input_format(const std::filesystem::path& p) {
  if (cooler::utils::is_cooler(p.string())) {
    return "cool";
  }
  if (cooler::utils::is_multires_file(p.string())) {
    return "mcool";
  }
  assert(hic::utils::is_hic_file(p));
  return "hic";
}

[[nodiscard]] static std::string infer_output_format(const std::filesystem::path& p) {
  const auto ext = p.extension();
  if (ext == ".hic") {
    return "hic";
  }
  if (ext == ".mcool") {
    return "mcool";
  }
  if (ext == ".cool") {
    return "cool";
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("unable to infer output file format from file name {}"), p));
}

[[nodiscard]] static std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& p,
                                                                 std::string_view format) {
  if (format == "cool") {
    return {cooler::File::open_read_only(p.string()).bin_size()};
  }
  if (format == "mcool") {
    return cooler::utils::list_resolutions(p);
  }
  assert(format == "hic");
  return hic::utils::list_resolutions(p);
}

// NOLINTNEXTLINE(misc-no-recursion)
[[nodiscard]] static std::string infer_assembly(const std::filesystem::path& p,
                                                std::uint32_t resolution, std::string_view format) {
  if (format == "cool") {
    auto assembly = cooler::File::open_read_only(p.string()).attributes().assembly;
    return !!assembly ? *assembly : "unknown";
  }
  if (format == "mcool") {
    return infer_assembly(fmt::format(FMT_STRING("{}::/resolutions/{}"), p.string(), resolution),
                          resolution, "cool");
  }
  assert(format == "hic");
  return hic::HiCFile{p.string(), resolution}.assembly();
}

void Cli::transform_args_convert_subcommand() {
  auto& c = std::get<ConvertConfig>(this->_config);

  c.input_format = infer_input_format(c.path_to_input);
  if (c.output_format.empty()) {
    c.output_format = infer_output_format(c.path_to_output);
  }

  if (c.resolutions.empty()) {
    c.resolutions = list_resolutions(c.path_to_input, c.input_format);
  }

  c.normalization_methods = generate_norm_vect(c.normalization_methods_str);

  if (c.genome.empty()) {
    c.genome = infer_assembly(c.path_to_input, c.resolutions.back(), c.input_format);
  }

  if (c.norm_dset_names.empty()) {
    c.norm_dset_names.resize(c.normalization_methods.size());
    std::transform(c.normalization_methods.begin(), c.normalization_methods.end(),
                   c.norm_dset_names.begin(), [](const auto norm) { return fmt::to_string(norm); });
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);
  c.verbosity = static_cast<std::uint8_t>(spdlog::level::critical) - c.verbosity;

  if (c.tmp_dir.empty()) {
    c.tmp_dir = c.path_to_output.parent_path();
  }
}

void Cli::transform_args_dump_subcommand() {
  auto& c = std::get<DumpConfig>(this->_config);

  c.format = infer_input_format(c.uri);
  if (c.format == "hic" && c.resolution == 0) {
    assert(c.table == "chroms");
    c.resolution = hic::utils::list_resolutions(c.uri).back();
  }

  if (this->_cli.get_subcommand("dump")->get_option("--range2")->empty()) {
    c.range2 = c.range1;
  }
}

void Cli::transform_args() {
  switch (this->_subcommand) {
    case convert:
      this->transform_args_convert_subcommand();
      break;
    case dump:
      this->transform_args_dump_subcommand();
      break;
    case load:
      [[fallthrough]];
    case merge:
      [[fallthrough]];
    case validate:
      [[fallthrough]];
    case help:
      break;
  }
}

}  // namespace hictk::tools
