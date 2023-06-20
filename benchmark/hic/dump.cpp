// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <cstdint>

#include "hictk/hic.hpp"

using namespace hictk;

struct Config {
  std::string path_to_hic{};
  std::vector<std::uint32_t> resolutions{};

  std::size_t target_num_records{1'000'000};
  bool genome_wide{};
  std::vector<std::size_t> block_cache_sizes{500'000'000};
};

static void dump_genome_wide(const std::string& path_to_hic, std::uint32_t resolution,
                             std::size_t target_num_records, std::size_t block_cache_size) {
  hic::HiCFile hf(path_to_hic, resolution, hic::MatrixType::observed, hic::MatrixUnit::BP,
                  block_cache_size);
  auto sel = hf.fetch();

  const auto t0 = std::chrono::steady_clock::now();
  auto first = sel.begin<float>();
  auto last = sel.end<float>();

  std::size_t i = 0;
  for (; i < target_num_records && first != last; ++i) {
    std::ignore = ++first;
  }

  const auto t1 = std::chrono::steady_clock::now();

  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
      1.0e6;
  fmt::print(FMT_STRING("{}\t{}\t{}\t{}\t{}\n"), path_to_hic, resolution, i, block_cache_size,
             delta);
}

static void dump(const std::string& path_to_hic, std::uint32_t resolution,
                 std::size_t target_num_records, std::size_t block_cache_size) {
  hic::HiCFile hf(path_to_hic, resolution, hic::MatrixType::observed, hic::MatrixUnit::BP,
                  block_cache_size);

  auto sel = hf.fetch(hf.chromosomes().longest_chromosome().name());

  const auto t0 = std::chrono::steady_clock::now();
  auto first = sel.begin<float>();
  auto last = sel.end<float>();

  std::size_t i = 0;
  for (; i < target_num_records && first != last; ++i) {
    std::ignore = ++first;
  }

  const auto t1 = std::chrono::steady_clock::now();

  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
      1.0e6;
  fmt::print(FMT_STRING("{}\t{}\t{}\t{}\t{}\n"), path_to_hic, resolution, i, block_cache_size,
             delta);
}

int main(int argc, char** argv) {
  CLI::App cli{};

  Config c{};

  cli.add_option("hic", c.path_to_hic, "Path to a .hic file.");
  cli.add_option("resolution", c.resolutions, "Resolution in bp.");
  cli.add_option("--target-num-records", c.target_num_records, "")->capture_default_str();
  cli.add_option("--block-cache-size", c.block_cache_sizes, "")->capture_default_str();
  cli.add_flag("--genome-wide", c.genome_wide, "")->capture_default_str();

  try {
    cli.parse(argc, argv);

    fmt::print(FMT_STRING("file\tresolution\tnum_records\tblock_cache_size\ttime\n"));
    for (const auto& resolution : c.resolutions) {
      for (const auto& block_size : c.block_cache_sizes) {
        c.genome_wide
            ? dump_genome_wide(c.path_to_hic, resolution, c.target_num_records, block_size)
            : dump(c.path_to_hic, resolution, c.target_num_records, block_size);
      }
    }
  } catch (const CLI::ParseError& e) {
    return cli.exit(e);
  }
}
