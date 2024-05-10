// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <hictk/cooler/cooler.hpp>
#include <hictk/transformers/coarsen.hpp>
#include <random>
#include <string_view>
#include <vector>

using namespace hictk;

struct Config {
  std::filesystem::path uri{};
  std::uint32_t factor{2};
  std::size_t iterations{1};
};

using PixelBuffer = std::vector<ThinPixel<std::uint32_t>>;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("uri", config.uri, "URI to a cooler file.")->required();
  cli.add_option("--factor", config.factor, "Zoomify factor.")->capture_default_str();
  cli.add_option("--iterations", config.iterations, "Number of iterations to perform.")
      ->capture_default_str();

  try {
    cli.parse(argc, argv);

    const cooler::File f(config.uri.string());
    const PixelBuffer pixels(f.fetch().begin<std::uint32_t>(), f.fetch().end<std::uint32_t>());
    PixelBuffer coarsened_pixels(pixels.size());
    coarsened_pixels.clear();

    std::uint64_t elapsed_time = 0;
    for (std::size_t i = 0; i < config.iterations; ++i) {
      const auto t0 = std::chrono::system_clock::now();
      const transformers::CoarsenPixels coarsener(pixels.begin(), pixels.end(), f.bins_ptr(),
                                                  config.factor);
      std::for_each(coarsener.begin(), coarsener.end(),
                    [&](const auto &p) { coarsened_pixels.push_back(p); });
      const auto t1 = std::chrono::system_clock::now();

      const auto delta = static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
      elapsed_time += delta;
    }

    const auto avg_time =
        (static_cast<double>(elapsed_time) / static_cast<double>(config.iterations)) / 1.0e9;
    const auto throughput = static_cast<double>(pixels.size()) / avg_time;

    fmt::print(FMT_STRING("hictk::transformers::Coarsener throughput: {:.4} pixels/s\n"),
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
