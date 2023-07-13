// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/cli.hpp"

#include <fmt/format.h>
#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <string>

namespace hictk::tools {

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
    } else if (this->_cli.get_subcommand("zoomify")->parsed()) {
      this->_subcommand = subcommand::zoomify;
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
    case zoomify:
      return "zoomify";
    default:
      assert(s == help);
      return "--help";
  }
}

void Cli::make_cli() {
  this->_cli.name(this->_exec_name);
  this->_cli.description("Coolerpp tools.");
  this->_cli.set_version_flag("-V,--version", std::string{hictk::config::version::str_long()});
  this->_cli.require_subcommand(1);

  this->make_convert_subcommand();
  this->make_dump_subcommand();
  this->make_load_subcommand();
  this->make_merge_subcommand();
  this->make_zoomify_subcommand();
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
    case zoomify:
      this->validate_zoomify_subcommand();
      break;
    case help:
      break;
  }
}

void Cli::transform_args() {
  switch (this->_subcommand) {
    case convert:
      this->transform_args_convert_subcommand();
      break;
    case dump:  // NOLINT
      this->transform_args_dump_subcommand();
      break;
    case load:   // NOLINT
      [[fallthrough]];
    case merge:  // NOLINT
      break;
    case zoomify:
      this->transform_args_zoomify_subcommand();
      break;
    case help:
      break;
  }
}

}  // namespace hictk::tools
