// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <hictk/bin.hpp>
#include <hictk/bin_table.hpp>
#include <hictk/chromosome.hpp>
#include <hictk/reference.hpp>
#include <random>
#include <string_view>
#include <vector>

using namespace hictk;

// clang-format off
const std::vector<Chromosome> hg38{
   Chromosome{0,  "chr1",  248956422},
   Chromosome{1,  "chr2",  242193529},
   Chromosome{2,  "chr3",  198295559},
   Chromosome{3,  "chr4",  190214555},
   Chromosome{4,  "chr5",  181538259},
   Chromosome{5,  "chr6",  170805979},
   Chromosome{6,  "chr7",  159345973},
   Chromosome{7,  "chr8",  145138636},
   Chromosome{8,  "chr9",  138394717},
   Chromosome{9,  "chr10", 133797422},
   Chromosome{10, "chr11", 135086622},
   Chromosome{11, "chr12", 133275309},
   Chromosome{12, "chr13", 114364328},
   Chromosome{13, "chr14", 107043718},
   Chromosome{14, "chr15", 101991189},
   Chromosome{15, "chr16", 90338345},
   Chromosome{16, "chr17", 83257441},
   Chromosome{17, "chr18", 80373285},
   Chromosome{18, "chr19", 58617616},
   Chromosome{19, "chr20", 64444167},
   Chromosome{20, "chr21", 46709983},
   Chromosome{21, "chr22", 50818468},
   Chromosome{22, "chrX",  156040895},
   Chromosome{23, "chrY",  57227415}
};

// clang-format on

struct Config {
  std::uint32_t resolution{1'000};
  std::size_t batch_size{10'000'000};
  std::size_t iterations{1};
  std::uint64_t seed{123456789};
};

[[nodiscard]] static std::vector<std::uint64_t> init_bin_ids(const BinTable &bins,
                                                             std::size_t batch_size,
                                                             std::uint64_t seed) {
  std::vector<std::uint64_t> buff(batch_size);
  std::mt19937_64 rand_eng(seed);
  std::generate(buff.begin(), buff.end(), [&]() {
    return std::uniform_int_distribution<std::uint64_t>{0, bins.size() - 1}(rand_eng);
  });

  return buff;
}

[[nodiscard]] static std::vector<Bin> init_bins(const BinTable &bins,
                                                const std::vector<std::uint64_t> &bin_ids) {
  std::vector<Bin> buff(bin_ids.size());
  for (std::size_t i = 0; i < bin_ids.size(); ++i) {
    buff[i] = bins.at(bin_ids[i]);
  }
  return buff;
}

[[nodiscard]] std::uint64_t run_benchmark(const BinTable &bins, const std::vector<Bin> &queries) {
  const auto t0 = std::chrono::system_clock::now();
  for (const auto &b : queries) {
    std::ignore = bins.at(b.chrom().name(), b.start());
  }
  const auto t1 = std::chrono::system_clock::now();

  return static_cast<std ::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("--resolution", config.resolution, "Resolution of the bin table.")
      ->capture_default_str();
  cli.add_option("--batch-size", config.batch_size, "Batch size.")->capture_default_str();
  cli.add_option("--iterations", config.iterations, "Number of iterations to perform.")
      ->capture_default_str();
  cli.add_option("--seed", config.seed, "Seed")->capture_default_str();

  try {
    cli.parse(argc, argv);

    const BinTable bin_table{hg38.begin(), hg38.end(), config.resolution};
    const auto bin_ids = init_bin_ids(bin_table, config.batch_size, config.seed);
    const auto bins = init_bins(bin_table, bin_ids);

    std::uint64_t elapsed_time = 0;
    for (std::size_t i = 0; i < config.iterations; ++i) {
      elapsed_time += run_benchmark(bin_table, bins);
    }

    const auto elapsed_seconds = static_cast<double>(elapsed_time) / 1.0e9;
    const auto throughput =
        static_cast<double>(config.batch_size * config.iterations) / elapsed_seconds;

    fmt::print(FMT_STRING("hictk::BinTable::at(chrom, pos) throughput: {:.4} num/s\n"), throughput);

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
