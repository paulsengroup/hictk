// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/cli.hpp"

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>

#include "hictk/tools/config.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}

auto Cli::parse_arguments() -> Config {
  try {
    _cli.name(_exec_name);
    _cli.parse(_argc, _argv);

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
  switch (s) {
    case balance:
      return "balance";
    case convert:
      return "convert";
    case dump:
      return "dump";
    case fix_mcool:
      return "fix-mcool";
    case load:
      return "load";
    case merge:
      return "merge";
    case metadata:
      return "metadata";
    case rename_chromosomes:
      return "rename-chromosomes";
    case validate:
      return "validate";
    case zoomify:
      return "zoomify";
    default:
      assert(s == help);
      return "--help";
  }
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description("Blazing fast tools to work with .hic and .cool files.");
  _cli.set_version_flag("-V,--version", std::string{hictk::config::version::str_long()});
  _cli.require_subcommand(1);

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
  switch (_subcommand) {
    case balance:
      validate_balance_subcommand();
      break;
    case convert:
      validate_convert_subcommand();
      break;
    case dump:
      validate_dump_subcommand();
      break;
    case fix_mcool:
      validate_fix_mcool_subcommand();
      break;
    case load:
      validate_load_subcommand();
      break;
    case merge:
      validate_merge_subcommand();
      break;
    case metadata:
      break;
    case rename_chromosomes:
      validate_rename_chromosomes_subcommand();
      break;
    case validate:
      break;
    case zoomify:
      validate_zoomify_subcommand();
      break;
    case help:
      break;
  }
}

void Cli::transform_args() {
  switch (_subcommand) {
    case balance:
      transform_args_balance_subcommand();
      break;
    case convert:
      transform_args_convert_subcommand();
      break;
    case dump:
      transform_args_dump_subcommand();
      break;
    case fix_mcool:
      transform_args_fix_mcool_subcommand();
      break;
    case load:
      transform_args_load_subcommand();
      break;
    case merge:
      transform_args_merge_subcommand();
      break;
    case metadata:
      transform_args_metadata_subcommand();
      break;
    case rename_chromosomes:
      transform_args_rename_chromosomes_subcommand();
      break;
    case validate:
      transform_args_validate_subcommand();
      break;
    case zoomify:
      transform_args_zoomify_subcommand();
      break;
    case help:
      break;
  }
}

}  // namespace hictk::tools
