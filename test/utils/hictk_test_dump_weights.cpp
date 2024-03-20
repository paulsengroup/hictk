// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <string>

#include "hictk/file.hpp"

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char** argv) noexcept {
  std::filesystem::path file{};
  std::uint32_t resolution{};
  std::string normalization{};
  bool divisive_weights{false};
  try {
    CLI::App cli{};

    cli.name(argv[0]);  // NOLINT

    cli.add_option("file", file, "Path to the input file.")->required();
    cli.add_option("--resolution", resolution, "Resolution in bp.")->required();
    cli.add_option("--normalization", normalization, "Normalization name.")->required();
    cli.add_flag("--divisive-weights", divisive_weights, "Return divisive weights.")
        ->capture_default_str();

    cli.parse(argc, argv);

    const hictk::File f(file, resolution);
    const auto weights = f.normalization(normalization);

    const auto weights_ = divisive_weights ? weights(hictk::balancing::Weights::Type::DIVISIVE)
                                           : weights(weights.type());

    fmt::print(FMT_STRING("{}\n"), fmt::join(weights_, "\n"));

  } catch (const std::exception& e) {
    fmt::print(stderr, FMT_STRING("error occurred while dumping \"{}\" weights from file {}: {}"),
               normalization, file, e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr, FMT_STRING("error occurred while dumping \"{}\" weights from file {}"),
               normalization, file);
    return 1;
  }
}
