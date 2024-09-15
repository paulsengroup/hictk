// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <cstdio>
#include <hictk/cooler/cooler.hpp>
#include <hictk/fmt/pixel.hpp>

using namespace hictk;

struct Config {
  std::filesystem::path uri{};
  bool join{false};
  std::size_t iterations{1};
};

template <typename PixelT>
[[nodiscard]] static std::ptrdiff_t print_pixels(const std::vector<PixelT> &pixels) {
  auto *dev_null = std::fopen("/dev/null", "w");
  std::for_each(pixels.begin(), pixels.end(),
                [&](const auto &p) { fmt::print(dev_null, FMT_COMPILE("{}\n"), p); });
  std::fclose(dev_null);

  return static_cast<std::ptrdiff_t>(pixels.size());
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("uri", config.uri, "URI to a cooler file.");
  cli.add_flag("--bg2,!--coo", config.join, "Join genomic coordinates.")->capture_default_str();
  cli.add_option("--iterations", config.iterations, "Number of iterations.")->capture_default_str();

  try {
    cli.parse(argc, argv);

    cooler::File f(config.uri.string());

    std::ptrdiff_t size = 0;
    std::uint64_t elapsed_time{};

    std::vector<Pixel<std::uint32_t>> pixel_buffer{};
    std::vector<ThinPixel<std::uint32_t>> thin_pixel_buffer{};

    std::for_each(f.begin<std::uint32_t>(), f.end<std::uint32_t>(), [&](const auto &tp) {
      if (config.join) {
        pixel_buffer.emplace_back(Pixel<std::uint32_t>(f.bins(), tp));

      } else {
        thin_pixel_buffer.emplace_back(tp);
      }
    });

    for (std::size_t i = 0; i < config.iterations; ++i) {
      const auto t0 = std::chrono::system_clock::now();
      size += print_pixels(thin_pixel_buffer);
      size += print_pixels(pixel_buffer);
      const auto t1 = std::chrono::system_clock::now();
      const auto delta = static_cast<std ::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
      elapsed_time += delta;
    }

    const auto elapsed_seconds = static_cast<double>(elapsed_time) / 1.0e9;
    const auto throughput = static_cast<double>(size) / elapsed_seconds;

    fmt::print(FMT_STRING("fmt::print({}) throughput: {:.4} num/s\n"),
               config.join ? "Pixel<std::uint32_t>" : "ThinPixel<std::uint32_t>", throughput);

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
