// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cassert>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string_view>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/tools.hpp"

using namespace hictk::tools;

int main(int argc, char** argv) {
  std::unique_ptr<Cli> cli{nullptr};
  try {
    cli = std::make_unique<Cli>(argc, argv);
    const auto config = cli->parse_arguments();

    {
      using subcmd = Cli::subcommand;
      switch (cli->get_subcommand()) {
        case subcmd::dump:
          dump_subcmd(std::get<DumpConfig>(config));
          break;
        case subcmd::load:
          // load_subcmd(std::get<LoadConfig>(config));
          break;
        case subcmd::merge:
          // merge_subcmd(std::get<MergeConfig>(config));
          break;
        default:
          throw std::runtime_error(
              "Default branch in switch statement in coolerpp_tools::main() should be unreachable! "
              "If "
              "you see this message, please file an issue on GitHub");
      }
    }
  } catch (const CLI::ParseError& e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    fmt::print(stderr, FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), err.what());
    return 1;
  } catch (const std::exception& e) {
    assert(cli);
    fmt::print(stderr,
               FMT_STRING("FAILURE! coolerpp_tools {} encountered the following error: {}."),
               cli->get_printable_subcommand(), e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! coolerpp_tools {} encountered the following error: Caught an "
                          "unhandled exception! "
                          "If you see this message, please file an issue on GitHub."),
               cli->get_printable_subcommand());
    return 1;
  }
  return 0;
}
