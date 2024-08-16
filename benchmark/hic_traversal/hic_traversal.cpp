// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <hictk/hic.hpp>

using namespace hictk;

struct Config {
  std::filesystem::path uri{};
  std::uint32_t resolution{};
  bool sorted{true};
  std::size_t iterations{1};
};

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("uri", config.uri, "Path to a .hic file.");
  cli.add_option("resolution", config.resolution, "Hi-C matrix resolution.");
  cli.add_option("--iterations", config.iterations, "Number of iterations.")->capture_default_str();
  cli.add_flag("--sorted,!--unsorted", config.sorted);

  try {
    cli.parse(argc, argv);

    hic::File f(config.uri.string(), config.resolution);

    std::ptrdiff_t size = 0;
    std::uint64_t elapsed_time{};

    for (std::size_t i = 0; i < config.iterations; ++i) {
      const auto t0 = std::chrono::system_clock::now();
      const auto sel = f.fetch();
      std::for_each(sel.begin<std::uint32_t>(config.sorted), sel.end<std::uint32_t>(),
                    [&]([[maybe_unused]] const auto n) { size++; });
      const auto t1 = std::chrono::system_clock::now();
      const auto delta = static_cast<std ::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
      elapsed_time += delta;
    }

    const auto elapsed_seconds = static_cast<double>(elapsed_time) / 1.0e9;
    const auto throughput = static_cast<double>(size) / elapsed_seconds;

    fmt::print(FMT_STRING("hictk::hic::File::iterator<std::uint32_t> throughput: {:.4} num/s\n"),
               throughput);

  } catch (const CLI::ParseError &e) {
    return cli.exit(e);
  } catch (const std::exception &e) {
    assert(cli);
    fmt::print(stderr, FMT_STRING("FAILURE! {} encountered the following error: {}.\n"), argv[0],
               e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! {} encountered the following error: Caught an "
                          "unhandled exception!\n"),
               argv[0]);
    return 1;
  }
  return 0;
}
