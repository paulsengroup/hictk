// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <hictk/cooler/cooler.hpp>
#include <hictk/hic/file_writer.hpp>
#include <hictk/pixel.hpp>
#include <hictk/tmpdir.hpp>
#include <vector>

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers)
struct Config {
  std::filesystem::path uri{};
  std::filesystem::path out_path{};
  std::size_t chunk_size{10'000'000};
  std::size_t threads{1};
  std::size_t iterations{1};
  bool validate{true};
};
// NOLINTEND(*-avoid-magic-numbers)

using PixelBuffer = std::vector<Pixel<float>>;

[[nodiscard]] static std::vector<PixelBuffer> chunk_pixels(const cooler::File &f,
                                                           const Reference &chroms,
                                                           std::size_t chunk_size) {
  std::vector<PixelBuffer> buffer{};
  PixelBuffer chunk{};

  const BinTable bins(chroms, f.resolution());

  std::for_each(f.begin<float>(), f.end<float>(), [&](const auto &p) {
    if (chunk.size() == chunk_size) {
      buffer.push_back(chunk);
      chunk.clear();
    }
    chunk.push_back(Pixel<float>(bins, p));
  });

  if (!chunk.empty()) {
    buffer.emplace_back(chunk);
  }

  return buffer;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  const auto *argv0 = argv[0];  // NOLINT(*-pointer-arithmetic)

  CLI::App cli{};
  Config config{};
  cli.add_option("in-uri", config.uri, "URI to an input cooler file.")->required();
  cli.add_option("out-path", config.out_path, "Path where to store the output cooler.")->required();
  cli.add_option("--chunk-size", config.chunk_size, "Chunk size.")->capture_default_str();
  cli.add_option("--iterations", config.iterations, "Number of iterations to perform.")
      ->capture_default_str();
  cli.add_option("--threads", config.threads, "Number of threads.");
  cli.add_flag("--validate,!--no-validate", config.validate, "Validate pixels before append.")
      ->capture_default_str();

  try {
    spdlog::default_logger()->set_level(spdlog::level::warn);
    cli.parse(argc, argv);

    const cooler::File f(config.uri.string());

    const auto chroms = f.chromosomes().remove_ALL();
    const auto pixels = chunk_pixels(f, chroms, config.chunk_size);

    std::size_t num_pixels = 0;
    std::uint64_t elapsed_time = 0;
    for (std::size_t i = 0; i < config.iterations; ++i) {
      const auto t0 = std::chrono::system_clock::now();
      {
        // NOLINTBEGIN(*-avoid-magic-numbers)
        hic::internal::HiCFileWriter writer(
            config.out_path.string(), chroms, {f.resolution()}, "unknown", config.threads,
            config.chunk_size, internal::TmpDir::default_temp_directory_path(), 11, true);
        // NOLINTEND(*-avoid-magic-numbers)
        for (const auto &chunk : pixels) {
          writer.add_pixels(f.resolution(), chunk.begin(), chunk.end());
          num_pixels += chunk.size();
        }
        writer.serialize();
      }
      const auto t1 = std::chrono::system_clock::now();
      std::filesystem::remove(config.out_path);  // NOLINT

      const auto delta = static_cast<std::uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
      elapsed_time += delta;
    }

    const auto avg_time =
        (static_cast<double>(elapsed_time) / static_cast<double>(config.iterations)) / 1.0e9;
    const auto throughput = static_cast<double>(num_pixels) / avg_time;

    fmt::print(FMT_STRING("hictk::hic::internal::HiCFileWriter throughput: {:.4} pixels/s\n"),
               throughput);

  } catch (const CLI::ParseError &e) {
    assert(cli);
    return cli.exit(e);
  } catch (const std::exception &e) {
    fmt::print(stderr, FMT_STRING("FAILURE! {} encountered the following error: {}.\n"), argv0,
               e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! {} encountered the following error: Caught an "
                          "unhandled exception!\n"),
               argv0);
    return 1;
  }
  return 0;
}
