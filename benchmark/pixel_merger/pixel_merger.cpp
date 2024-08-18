// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <hictk/pixel.hpp>
#include <hictk/transformers/pixel_merger.hpp>
#include <random>
#include <string_view>
#include <vector>

using namespace hictk;

struct Config {
  std::size_t genome_size{3'300'000};
  std::size_t num_pixels_per_chunk{100'000'000};
  std::size_t num_chunks{2};
  std::uint32_t resolution{1'000};
  std::size_t iterations{1};
  std::uint64_t seed{123456789};
};

using PixelBuffer = std::vector<ThinPixel<std::uint32_t>>;

[[nodiscard]] static PixelBuffer generate_pixels(std::size_t num_bins, std::size_t num_pixels,
                                                 std::mt19937_64 &rand_eng) {
  PixelBuffer buffer(num_pixels);

  std::generate(buffer.begin(), buffer.end(), [&]() {
    const auto bin1_id = std::uniform_int_distribution<std::uint64_t>{0, num_bins - 1}(rand_eng);
    const auto bin2_id =
        std::uniform_int_distribution<std::uint64_t>{bin1_id, num_bins - 1}(rand_eng);

    return ThinPixel<std::uint32_t>{bin1_id, bin2_id, 1};
  });

  std::sort(buffer.begin(), buffer.end());
  return buffer;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("--genome-size", config.genome_size, "Genome size in bp.")->capture_default_str();
  cli.add_option("--resolution", config.resolution, "Resolution of the bin table.")
      ->capture_default_str();
  cli.add_option("--pixels-per-chunk", config.num_pixels_per_chunk,
                 "Number of pixels found in each chunk to be merged.")
      ->capture_default_str();
  cli.add_option("--num-chunks", config.num_chunks, "Number of chunks to be merged.")
      ->capture_default_str();
  cli.add_option("--iterations", config.iterations, "Number of iterations to perform.")
      ->capture_default_str();
  cli.add_option("--seed", config.seed, "Seed")->capture_default_str();

  try {
    cli.parse(argc, argv);
    std::mt19937_64 rand_eng(config.seed);
    std::vector<PixelBuffer> pixel_chunks{};
    using PixelIt = decltype(pixel_chunks.front().begin());

    std::vector<PixelIt> heads{};
    std::vector<PixelIt> tails{};

    for (std::size_t i = 0; i < config.num_chunks; ++i) {
      pixel_chunks.emplace_back(generate_pixels(config.genome_size / config.resolution,
                                                config.num_pixels_per_chunk, rand_eng));
      heads.emplace_back(pixel_chunks.back().begin());
      tails.emplace_back(pixel_chunks.back().end());
    }

    std::uint64_t elapsed_time = 0;
    std::uint64_t sum = 0;
    for (std::size_t i = 0; i < config.iterations; ++i) {
      const transformers::PixelMerger merger(heads, tails);

      const auto t0 = std::chrono::system_clock::now();
      std::for_each(merger.begin(), merger.end(), [&](const auto &p) { sum += p.count; });
      const auto t1 = std::chrono::system_clock::now();

      const auto delta = static_cast<std ::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
      elapsed_time += delta;
    }

    std::ignore = sum;

    const auto avg_time =
        (static_cast<double>(elapsed_time) / static_cast<double>(config.iterations)) / 1.0e9;
    const auto throughput =
        static_cast<double>(config.num_chunks * config.num_pixels_per_chunk) / avg_time;

    fmt::print(FMT_STRING("hictk::transformers::PixelMerger throughput: {:.4} pixels/s\n"),
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
