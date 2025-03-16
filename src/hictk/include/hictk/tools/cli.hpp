// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <filesystem>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "config.hpp"
#include "hictk/cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/string_utils.hpp"

namespace hictk::hic {

inline bool lexical_cast(const std::string& input, MatrixType& v) {
  if (input == "observed") {
    v = MatrixType::observed;
    return true;
  }
  if (input == "oe") {
    v = MatrixType::oe;
    return true;
  }

  if (input == "expected") {
    v = MatrixType::expected;
    return true;
  }

  std::string input_lower{input};
  std::transform(input_lower.begin(), input_lower.end(), input_lower.begin(),
                 [](auto c) { return std::tolower(c); });

  if (input_lower != input) {
    return lexical_cast(input_lower, v);
  }

  return false;
}

inline bool lexical_cast(const std::string& input, MatrixUnit& v) {
  if (input == "BP") {
    v = MatrixUnit::BP;
    return true;
  }
  if (input == "FRAG") {
    v = MatrixUnit::FRAG;
    return true;
  }

  std::string input_upper{input};
  std::transform(input_upper.begin(), input_upper.end(), input_upper.begin(),
                 [](auto c) { return std::toupper(c); });

  if (input_upper != input) {
    return lexical_cast(input_upper, v);
  }

  return false;
}

}  // namespace hictk::hic

namespace hictk::tools {

class CoolerFileValidator : public CLI::Validator {
 public:
  CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_cooler(uri)) {
        if (cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        if (cooler::utils::is_scool_file(uri)) {
          return "URI points to a .scool file: " + uri;
        }
        const auto path = cooler::parse_cooler_uri(uri).file_path;
        if (!std::filesystem::exists(path)) {
          return "No such file: " + path;
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
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!cooler::utils::is_multires_file(uri)) {
        return "Not a valid multi-resolution cooler: " + uri;
      }
      return "";
    };
  }
};

class SingleCellCoolerFileValidator : public CLI::Validator {
 public:
  SingleCellCoolerFileValidator() : Validator("Single-cell-cooler") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!cooler::utils::is_scool_file(uri)) {
        return "Not a valid single-cell cooler: " + uri;
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
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hic::utils::is_hic_file(path)) {
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
  [[nodiscard]] std::string make_option_opts(const CLI::Option* opt) const override {
    if (!opt->get_option_text().empty()) {
      return opt->get_option_text();
    }

    auto str_contains = [](const auto& s, const auto query) {
      return s.find(query) != remove_cvref_t<decltype(s)>::npos;
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

template <typename Enum>
class StringToEnumChecked : public CLI::Validator {
  static_assert(std::is_enum_v<Enum>);
  [[nodiscard]] static std::string to_lower(std::string_view s) {
    std::string s_lower{s};
    std::transform(s_lower.begin(), s_lower.end(), s_lower.begin(),
                   [](auto c) { return std::tolower(c); });
    return s_lower;
  }

 public:
  // NOLINTNEXTLINE(*-unnecessary-value-param)
  explicit StringToEnumChecked(std::vector<std::pair<std::string, Enum>> mappings) {
    assert(!mappings.empty());

    auto description_formatter = [mappings]() {
      std::vector<std::string> keys(mappings.size());
      std::transform(mappings.begin(), mappings.end(), keys.begin(),
                     [](const auto& kv) { return kv.first; });
      return fmt::format(FMT_STRING("{{{}}}"), fmt::join(keys, ","));
    };

    desc_function_ = description_formatter;

    func_ = [mappings, description_formatter](std::string& input) -> std::string {
      const auto input_lower = to_lower(input);

      const auto match = std::find_if(mappings.begin(), mappings.end(), [&](const auto& kv) {
        return to_lower(kv.first) == input_lower;
      });

      if (match != mappings.end()) {
        return "";
      }

      return fmt::format(FMT_STRING("{} not in {}"), input, description_formatter());
    };
  }
};

// clang-format off
inline const auto IsValidCoolerFile = CoolerFileValidator();                      // NOLINT(cert-err58-cpp)
inline const auto IsValidMultiresCoolerFile = MultiresCoolerFileValidator();      // NOLINT(cert-err58-cpp)
inline const auto IsValidSingleCellCoolerFile = SingleCellCoolerFileValidator();  // NOLINT(cert-err58-cpp)
inline const auto IsValidHiCFile = HiCFileValidator();                            // NOLINT(cert-err58-cpp)
// clang-format on

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
inline const auto ParseHiCMatrixType =
  StringToEnumChecked{
    std::vector<std::pair<std::string, hic::MatrixType>>{
      {"observed", hic::MatrixType::observed},
      {"oe", hic::MatrixType::oe},
      {"expected", hic::MatrixType::expected}
    }
  }.description("");

// NOLINTNEXTLINE(cert-err58-cpp)
inline const auto ParseHiCMatrixUnit =
  StringToEnumChecked{
    std::vector<std::pair<std::string, hic::MatrixUnit>>{
      {"BP", hic::MatrixUnit::BP},
      {"FRAG", hic::MatrixUnit::FRAG}
    }
  }.description("");
// clang-format on

class Cli {
 public:
  enum class subcommand : std::uint_fast8_t {
    help,
    balance,
    convert,
    dump,
    fix_mcool,
    load,
    merge,
    metadata,
    rename_chromosomes,
    validate,
    zoomify,
  };

  Cli(int argc, char** argv);
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] auto parse_arguments() -> Config;
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] int exit() const noexcept;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;
  void log_warnings() const noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  Config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};
  mutable std::vector<std::string> _warnings{};

  void make_balance_subcommand();
  void make_ice_balance_subcommand(CLI::App& app);
  void make_scale_balance_subcommand(CLI::App& app);
  void make_vc_balance_subcommand(CLI::App& app);
  void make_convert_subcommand();
  void make_dump_subcommand();
  void make_fix_mcool_subcommand();
  void make_load_subcommand();
  void make_merge_subcommand();
  void make_metadata_subcommand();
  void make_rename_chromosomes_subcommand();
  void make_validate_subcommand();
  void make_zoomify_subcommand();
  void make_cli();

  void validate_balance_subcommand() const;
  void validate_convert_subcommand() const;
  void validate_dump_subcommand() const;
  void validate_fix_mcool_subcommand() const;
  void validate_load_subcommand() const;
  void validate_merge_subcommand() const;
  void validate_rename_chromosomes_subcommand() const;
  void validate_zoomify_subcommand() const;
  void validate_args() const;

  void transform_args_balance_subcommand();
  void transform_args_ice_balance_subcommand();
  void transform_args_scale_balance_subcommand();
  void transform_args_vc_balance_subcommand();
  void transform_args_convert_subcommand();
  void transform_args_dump_subcommand();
  void transform_args_fix_mcool_subcommand();
  void transform_args_load_subcommand();
  void transform_args_merge_subcommand();
  void transform_args_metadata_subcommand();
  void transform_args_rename_chromosomes_subcommand();
  void transform_args_validate_subcommand();
  void transform_args_zoomify_subcommand();
  void transform_args();
};

[[nodiscard]] inline std::string infer_input_format(const std::filesystem::path& p) {
  if (hic::utils::is_hic_file(p)) {
    return "hic";
  }
  if (cooler::utils::is_cooler(p.string())) {
    return "cool";
  }
  if (cooler::utils::is_multires_file(p.string())) {
    return "mcool";
  }
  if (cooler::utils::is_scool_file(p.string())) {
    return "scool";
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("unable to infer file format for file \"{}\""), p.string()));
}

[[nodiscard]] inline std::string infer_output_format(const std::filesystem::path& p) {
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

[[nodiscard]] inline std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& p,
                                                                 std::string_view format) {
  if (format == "cool") {
    return {cooler::File(p.string()).resolution()};
  }
  if (format == "scool") {
    return {cooler::SingleCellFile{p.string()}.resolution()};
  }
  if (format == "mcool") {
    return cooler::utils::list_resolutions(p, true);
  }
  assert(format == "hic");
  return hic::utils::list_resolutions(p, true);
}

[[nodiscard]] inline std::optional<std::int16_t> parse_hictk_verbosity_from_env(
    bool skip = false) noexcept {
  if (skip) {
    return {};
  }

  constexpr auto critical = static_cast<std::int16_t>(spdlog::level::critical);

  if (std::getenv("HICTK_QUIET")) {  // NOLINT(*-mt-unsafe)
    return critical;
  }

  if (std::getenv("VERBOSE")) {  // NOLINT(*-mt-unsafe)
    return static_cast<std::int16_t>(spdlog::level::debug);
  }

  const auto* verbosity_ptr = std::getenv("HICTK_VERBOSITY");  // NOLINT(*-mt-unsafe)
  if (!verbosity_ptr) {
    return {};
  }

  std::string_view verbosity{verbosity_ptr};

  // in spdlog, high numbers correspond to low log levels
  if (verbosity == "0") {
    return critical;
  }
  if (verbosity == "1") {
    return critical - 1;
  }
  if (verbosity == "2") {
    return critical - 2;
  }
  if (verbosity == "3") {
    return critical - 3;
  }
  if (verbosity == "4") {
    return critical - 4;
  }
  if (verbosity == "5") {
    return critical - 5;
  }

  auto str_compare_case_insensitive = [&](std::string_view lvl) noexcept {
    if (verbosity.size() != lvl.size()) {
      return false;
    }

    for (std::size_t i = 0; i < verbosity.size(); ++i) {
      if (std::tolower(verbosity[i]) != lvl[i]) {
        return false;
      }
    }
    return true;
  };

  if (str_compare_case_insensitive("critical")) {
    return critical;
  }
  if (str_compare_case_insensitive("error") || str_compare_case_insensitive("err")) {
    return static_cast<std::int16_t>(spdlog::level::err);
  }
  if (str_compare_case_insensitive("warning") || str_compare_case_insensitive("warn")) {
    return static_cast<std::int16_t>(spdlog::level::warn);
  }
  if (str_compare_case_insensitive("info")) {
    return static_cast<std::int16_t>(spdlog::level::info);
  }
  if (str_compare_case_insensitive("debug")) {
    return static_cast<std::int16_t>(spdlog::level::debug);
  }

  fmt::println(
      stderr,
      FMT_STRING(
          "WARNING: unable to parse verbosity level from env variable HICTK_VERBOSITY=\"{}\""),
      verbosity);
  return {};
}

}  // namespace hictk::tools
