// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <CLI/CLI.hpp>
#include <string>
#include <string_view>

#include "config.hpp"

namespace hictk::tools {

class Cli {
 public:
  enum subcommand {
    help,
    convert,
    dump,
    load,
    merge,
  };
  Cli(int argc, char** argv);
  [[nodiscard]] subcommand get_subcommand() const noexcept;
  [[nodiscard]] std::string_view get_printable_subcommand() const noexcept;
  [[nodiscard]] auto parse_arguments() -> Config;
  [[nodiscard]] int exit(const CLI::ParseError& e) const;
  [[nodiscard]] static std::string_view subcommand_to_str(subcommand s) noexcept;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  Config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  void make_convert_subcommand();
  void make_dump_subcommand();
  void make_load_subcommand();
  void make_merge_subcommand();
  void make_cli();

  void validate_convert_subcommand() const;
  void validate_dump_subcommand() const;
  void validate_load_subcommand() const;
  void validate_merge_subcommand() const;
  void validate() const;

  void transform_args_convert_subcommand();
  void transform_args_dump_subcommand();
  void transform_args();
};

}  // namespace hictk::tools
