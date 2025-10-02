// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/config.hpp"

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

class Cli {
 public:
  enum class subcommand : std::uint_fast8_t {
    none,
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
  subcommand _subcommand{subcommand::none};
  mutable std::vector<std::string> _warnings{};
  std::string_view _help_flag{};

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

  [[nodiscard]] bool handle_help_flags();
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
  try {
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

    // NOLINTBEGIN(*-avoid-magic-numbers)
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
    // NOLINTEND(*-avoid-magic-numbers)

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
  } catch (...) {
    return {};
  }
}

}  // namespace hictk::tools
