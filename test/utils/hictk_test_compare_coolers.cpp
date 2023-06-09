// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <CLI/CLI.hpp>
#include <filesystem>

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char** argv) noexcept {
  std::filesystem::path file1{};
  std::filesystem::path file2{};
  try {
    CLI::App cli{};

    cli.name(argv[0]);  // NOLINT

    cli.add_option("cooler1", file1, "Path to the first cooler file to compare.")->required();
    cli.add_option("cooler2", file2, "Path to the second cooler file to compare.")->required();

    cli.parse(argc, argv);

    const auto equal = hictk::cooler::utils::equal(file1.string(), file2.string());

    if (equal) {
      fmt::print(stderr, FMT_STRING("files are equal!\n"));
      return 0;
    }

    fmt::print(stderr, FMT_STRING("files are different!\n"));
    return 1;
  } catch (const std::exception& e) {
    fmt::print(stderr, FMT_STRING("error occurred while comparing {} with {}: {}"), file1, file2,
               e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr, FMT_STRING("error occurred while comparing {} with {}"), file1, file2);
    return 1;
  }
}
