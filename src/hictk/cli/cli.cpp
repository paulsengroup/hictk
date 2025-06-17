// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/cli.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>

#include "hictk/license.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/telemetry.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}

void Cli::log_warnings() const noexcept {
  for (const auto& w : _warnings) {
    SPDLOG_WARN(FMT_STRING("{}"), w);
  }
  _warnings.clear();
}

auto Cli::parse_arguments() -> Config {
  try {
    _cli.name(_exec_name);
    _cli.parse(_argc, _argv);

    if (handle_help_flags()) {
      return _config;
    }

    if (_cli.get_subcommand("balance")->parsed()) {
      _subcommand = subcommand::balance;
    } else if (_cli.get_subcommand("convert")->parsed()) {
      _subcommand = subcommand::convert;
    } else if (_cli.get_subcommand("dump")->parsed()) {
      _subcommand = subcommand::dump;
    } else if (_cli.get_subcommand("fix-mcool")->parsed()) {
      _subcommand = subcommand::fix_mcool;
    } else if (_cli.get_subcommand("load")->parsed()) {
      _subcommand = subcommand::load;
    } else if (_cli.get_subcommand("merge")->parsed()) {
      _subcommand = subcommand::merge;
    } else if (_cli.get_subcommand("metadata")->parsed()) {
      _subcommand = subcommand::metadata;
    } else if (_cli.get_subcommand("rename-chromosomes")->parsed()) {
      _subcommand = subcommand::rename_chromosomes;
    } else if (_cli.get_subcommand("validate")->parsed()) {
      _subcommand = subcommand::validate;
    } else if (_cli.get_subcommand("zoomify")->parsed()) {
      _subcommand = subcommand::zoomify;
    } else {
      _subcommand = subcommand::none;
      for (const auto* opt : {"--help", "--version"}) {
        if (!_cli.get_option(opt)->empty()) {
          _exit_code = 0;
          return _config;
        }
      }
      fmt::print(stderr, FMT_STRING("A subcommand is required\n"
                                    "Run with --help for more information.\n"));
      _exit_code = 1;
      return _config;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    _exit_code = exit(e);
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

int Cli::exit(const CLI::ParseError& e) const { return _cli.exit(e) != 0; }  // NOLINT
int Cli::exit() const noexcept { return _exit_code; }

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  using sc = subcommand;
  switch (s) {
    case sc::balance:
      return "balance";
    case sc::convert:
      return "convert";
    case sc::dump:
      return "dump";
    case sc::fix_mcool:
      return "fix-mcool";
    case sc::load:
      return "load";
    case sc::merge:
      return "merge";
    case sc::metadata:
      return "metadata";
    case sc::rename_chromosomes:
      return "rename-chromosomes";
    case sc::validate:
      return "validate";
    case sc::zoomify:
      return "zoomify";
    default:
      assert(s == sc::none);
      return "";
  }
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description("Blazing fast tools to work with .hic and .cool files.");
  _cli.set_version_flag("-V,--version", std::string{config::version::str_long()});

  auto* grp = _cli.add_option_group("help");
  grp->require_option(0, 1);

  grp->add_flag_callback(
      "--help-cite", [this]() { _help_flag = "cite"; },
      "Print hictk's citation in Bibtex format and exit.");
  grp->add_flag_callback(
      "--help-docs", [this]() { _help_flag = "docs"; },
      "Print the URL to hictk's documentation and exit.");
  grp->add_flag_callback(
      "--help-license", [this]() { _help_flag = "license"; }, "Print the hictk license and exit.");
  grp->add_flag_callback(
      "--help-telemetry", [this]() { _help_flag = "telemetry"; },
      "Print information regarding telemetry collection and exit.");

  make_balance_subcommand();
  make_convert_subcommand();
  make_dump_subcommand();
  make_fix_mcool_subcommand();
  make_load_subcommand();
  make_merge_subcommand();
  make_metadata_subcommand();
  make_rename_chromosomes_subcommand();
  make_validate_subcommand();
  make_zoomify_subcommand();
}

void Cli::validate_args() const {
  using sc = subcommand;
  switch (_subcommand) {
    case sc::balance:
      validate_balance_subcommand();
      break;
    case sc::convert:
      validate_convert_subcommand();
      break;
    case sc::dump:
      validate_dump_subcommand();
      break;
    case sc::fix_mcool:
      validate_fix_mcool_subcommand();
      break;
    case sc::load:
      validate_load_subcommand();
      break;
    case sc::merge:
      validate_merge_subcommand();
      break;
    case sc::metadata:
      break;
    case sc::rename_chromosomes:
      validate_rename_chromosomes_subcommand();
      break;
    case sc::validate:
      break;
    case sc::zoomify:
      validate_zoomify_subcommand();
      break;
    case sc::none:
      break;
  }
}

void Cli::transform_args() {
  using sc = subcommand;
  switch (_subcommand) {
    case sc::balance:
      transform_args_balance_subcommand();
      break;
    case sc::convert:
      transform_args_convert_subcommand();
      break;
    case sc::dump:
      transform_args_dump_subcommand();
      break;
    case sc::fix_mcool:
      transform_args_fix_mcool_subcommand();
      break;
    case sc::load:
      transform_args_load_subcommand();
      break;
    case sc::merge:
      transform_args_merge_subcommand();
      break;
    case sc::metadata:
      transform_args_metadata_subcommand();
      break;
    case sc::rename_chromosomes:
      transform_args_rename_chromosomes_subcommand();
      break;
    case sc::validate:
      transform_args_validate_subcommand();
      break;
    case sc::zoomify:
      transform_args_zoomify_subcommand();
      break;
    case sc::none:
      break;
  }
}

[[nodiscard]] static constexpr std::string_view get_citation() noexcept {
  // clang-format off
  return
R"(@article{hictk,
    author = {Rossini, Roberto and Paulsen, Jonas},
    title = "{hictk: blazing fast toolkit to work with .hic and .cool files}",
    journal = {Bioinformatics},
    volume = {40},
    number = {7},
    pages = {btae408},
    year = {2024},
    month = {06},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btae408},
    url = {https://doi.org/10.1093/bioinformatics/btae408},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/7/btae408/58385157/btae408.pdf},
}
)";
  // clang-format on
}

[[nodiscard]] static std::string_view get_telemetry_help() noexcept {
  // clang-format off
  if constexpr (Tracer::is_dummy()) {
    assert(!Tracer::should_collect_telemetry());
    return
  R"(hictk was compiled WITHOUT support for telemetry.
No telemetry data will be collected.
See https://hictk.readthedocs.io/en/stable/telemetry.html for more details.
)";
  }

  if (Tracer::should_collect_telemetry()) {
  return
R"(hictk was compiled WITH support for telemetry.
Telemetry data will be collected as the environment variable "HICTK_NO_TELEMETRY" is not defined.
See https://hictk.readthedocs.io/en/stable/telemetry.html for more details.
)";
  }
  return
R"(hictk was compiled WITH support for telemetry.
Telemetry data won't be collected as the environment variable "HICTK_NO_TELEMETRY" is defined.
See https://hictk.readthedocs.io/en/stable/telemetry.html for more details.
)";
  // clang-format on
}

bool Cli::handle_help_flags() {
  if (_help_flag.empty()) {
    return false;
  }

  if (_help_flag == "telemetry") {
    fmt::print(FMT_STRING("{}"), get_telemetry_help());
    _subcommand = subcommand::none;
    _exit_code = 0;
    return true;
  }

  if (_help_flag == "license") {
    fmt::print(FMT_STRING("{}"), config::license::license);
    _subcommand = subcommand::none;
    _exit_code = 0;
    return true;
  }

  if (_help_flag == "cite") {
    fmt::print(FMT_STRING("{}"), get_citation());
    _subcommand = subcommand::none;
    _exit_code = 0;
    return true;
  }

  if (_help_flag == "docs") {
    fmt::println(FMT_STRING("https://hictk.readthedocs.io"));
    _subcommand = subcommand::none;
    _exit_code = 0;
    return true;
  }

  HICTK_UNREACHABLE_CODE;
  return false;
}

}  // namespace hictk::tools
