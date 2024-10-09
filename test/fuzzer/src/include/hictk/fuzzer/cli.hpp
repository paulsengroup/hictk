// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <CLI/CLI.hpp>
#include <cassert>
#include <filesystem>
#include <map>
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

namespace hictk::fuzzer {

class CoolerFileValidator : public CLI::Validator {
 public:
  inline CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_cooler(uri)) {
        if (hictk::cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        if (hictk::cooler::utils::is_scool_file(uri)) {
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
  inline MultiresCoolerFileValidator() : Validator("Multires-cooler") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hictk::cooler::utils::is_multires_file(uri)) {
        return "Not a valid multi-resolution cooler: " + uri;
      }
      return "";
    };
  }
};

class HiCFileValidator : public CLI::Validator {
 public:
  inline HiCFileValidator() : Validator("HiC") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
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
inline const auto IsValidCoolerFile = CoolerFileValidator();                      // NOLINT(cert-err58-cpp)
inline const auto IsValidMultiresCoolerFile = MultiresCoolerFileValidator();      // NOLINT(cert-err58-cpp)
inline const auto IsValidHiCFile = HiCFileValidator();                            // NOLINT(cert-err58-cpp)
// clang-format on

class Cli {
 public:
  enum class subcommand : std::uint_fast8_t {
    help,
    fuzz,
    launch_worker,
  };

  Cli(int argc, char** argv);
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] auto parse_arguments() -> Config;
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] int exit() const noexcept;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  Config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  void make_fuzz_subcommand();
  void make_launch_worker_subcommand();
  void make_cli();

  void validate_fuzz_subcommand() const;
  void validate_launch_worker_subcommand() const;
  void validate_args() const;

  void transform_args_fuzz_subcommand();
  void transform_args_launch_worker_subcommand();
  void transform_args();
};

}  // namespace hictk::fuzzer
