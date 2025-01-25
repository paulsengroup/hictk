// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/benchmark/catch_constructor.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "./common.hpp"
#include "hictk/bin.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/bin_table_fixed.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

using namespace hictk;

using N = std::uint32_t;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
[[nodiscard]] static std::pair<Bin, Bin> sample_bin_pair(const BinTableFixed& bin_table,
                                                         std::mt19937_64& rand_eng, bool cis_pair) {
  const auto chrom1 = bin_table.chromosomes().at(std::uniform_int_distribution<std::uint32_t>{
      0, static_cast<std::uint32_t>(bin_table.num_chromosomes() - 1)}(rand_eng));

  auto sample_trans_chrom = [&]() {
    while (true) {
      auto chrom = bin_table.chromosomes().at(std::uniform_int_distribution<std::uint32_t>{
          0, static_cast<std::uint32_t>(bin_table.num_chromosomes() - 1)}(rand_eng));
      if (chrom != chrom1) {
        return chrom;
      }
    }
  };

  const auto chrom2 = cis_pair ? chrom1 : sample_trans_chrom();

  auto bin1 = bin_table.at(
      chrom1, std::uniform_int_distribution<std::uint32_t>{0, chrom1.size() - 1}(rand_eng));
  auto bin2 = bin_table.at(
      chrom2, std::uniform_int_distribution<std::uint32_t>{0, chrom2.size() - 1}(rand_eng));
  if (bin1 > bin2) {
    std::swap(bin1, bin2);
  }

  return std::make_pair(bin1, bin2);
}

[[nodiscard]] static std::vector<ThinPixel<N>> generate_cis_thin_pixels(
    const BinTableFixed& bin_table, std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<ThinPixel<N>> pixels(size);

  std::generate(pixels.begin(), pixels.end(), [&]() {
    const auto [bin1, bin2] = sample_bin_pair(bin_table, rand_eng, true);
    const auto count = std::uniform_int_distribution<std::uint32_t>{1, 1'000'000}(rand_eng);
    return ThinPixel<N>{bin1.id(), bin2.id(), count};
  });

  return pixels;
}

[[nodiscard]] static std::vector<ThinPixel<N>> generate_trans_thin_pixels(
    const BinTableFixed& bin_table, std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<ThinPixel<N>> pixels(size);

  std::generate(pixels.begin(), pixels.end(), [&]() {
    const auto [bin1, bin2] = sample_bin_pair(bin_table, rand_eng, false);
    const auto count = std::uniform_int_distribution<std::uint32_t>{1, 1'000'000}(rand_eng);
    return ThinPixel<N>{bin1.id(), bin2.id(), count};
  });

  return pixels;
}

[[nodiscard]] static std::vector<Pixel<N>> generate_cis_pixels(const BinTableFixed& bin_table,
                                                               std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<Pixel<N>> pixels(size);

  std::generate(pixels.begin(), pixels.end(), [&]() {
    const auto [bin1, bin2] = sample_bin_pair(bin_table, rand_eng, true);
    const auto count = std::uniform_int_distribution<std::uint32_t>{1, 1'000'000}(rand_eng);
    return Pixel{bin1, bin2, count};
  });

  return pixels;
}

[[nodiscard]] static std::vector<Pixel<N>> generate_trans_pixels(const BinTableFixed& bin_table,
                                                                 std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<Pixel<N>> pixels(size);

  std::generate(pixels.begin(), pixels.end(), [&]() {
    const auto [bin1, bin2] = sample_bin_pair(bin_table, rand_eng, false);
    const auto count = std::uniform_int_distribution<std::uint32_t>{1, 1'000'000}(rand_eng);
    return Pixel{bin1, bin2, count};
  });

  return pixels;
}

[[nodiscard]] static std::vector<ThinPixel<N>> generate_thin_pixels(std::size_t size) {
  const BinTableFixed bin_table({hg38.begin(), hg38.end()}, 1'000);
  const auto num_cis_interactions = static_cast<std::size_t>(0.7 * static_cast<double>(size));

  auto cis_pixels = generate_cis_thin_pixels(bin_table, num_cis_interactions);
  const auto trans_pixels = generate_trans_thin_pixels(bin_table, size - num_cis_interactions);

  cis_pixels.insert(cis_pixels.end(), trans_pixels.begin(), trans_pixels.end());
  return cis_pixels;
}

[[nodiscard]] static std::vector<Pixel<N>> generate_pixels(std::size_t size) {
  const BinTableFixed bin_table({hg38.begin(), hg38.end()}, 1'000);
  const auto num_cis_interactions = static_cast<std::size_t>(0.7 * static_cast<double>(size));

  auto cis_pixels = generate_cis_pixels(bin_table, num_cis_interactions);
  const auto trans_pixels = generate_trans_pixels(bin_table, size - num_cis_interactions);

  cis_pixels.insert(cis_pixels.end(), trans_pixels.begin(), trans_pixels.end());
  return cis_pixels;
}

TEST_CASE("Pixel") {
  BENCHMARK_ADVANCED("Construction")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = static_cast<std::size_t>(meter.runs());
    std::vector<Catch::Benchmark::storage_for<Pixel<N>>> storage(num_runs);

    const Chromosome chrom{0, "chr1", 123'456'789};
    const Bin bin{0, 0, chrom, 0, 1'000};

    meter.measure([&storage, &bin](std::size_t i) { storage[i].construct(bin, bin, 1); });
  };

  BENCHMARK_ADVANCED("Destruction")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = static_cast<std::size_t>(meter.runs());
    std::vector<Catch::Benchmark::destructable_object<Pixel<N>>> storage(num_runs);

    const Chromosome chrom{0, "chr1", 123'456'789};
    const Bin bin{0, 0, chrom, 0, 1'000};

    for (auto& pixel : storage) {
      pixel.construct(bin, bin, 1);
    }

    meter.measure([&storage](std::size_t i) { storage[i].destruct(); });
  };

  BENCHMARK_ADVANCED("from_coo (uint32)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<std::uint32_t>::from_coo(bin_table, "123456\t234567\t123");
    });
  };

  BENCHMARK_ADVANCED("from_coo (double)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure(
        [&bin_table]() { return Pixel<double>::from_coo(bin_table, "123456\t234567\t123.4567"); });
  };

  BENCHMARK_ADVANCED("from_bg2 (uint32)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<std::uint32_t>::from_bg2(bin_table,
                                            "chr7\t1000000\t1001000\tchr12\t1000000\t1001000\t123");
    });
  };

  BENCHMARK_ADVANCED("from_bg2 (double)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<double>::from_bg2(bin_table,
                                     "chr7\t1000000\t1001000\tchr12\t1000000\t1001000\t123.4567");
    });
  };

  BENCHMARK_ADVANCED("from_validpair")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<std::uint32_t>::from_validpair(
          bin_table,
          "NS500537:79:HFYYWBGX2:1:11112:2304:13920\tchr2\t12233\t+\tchr2\t13674\t+"
          "\t1\tfrag1\tfrag2\t1\t1\tallele-info");
    });
  };

  BENCHMARK_ADVANCED("from_4dn_pairs")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<std::uint32_t>::from_4dn_pairs(
          bin_table,
          "NS500537:79:HFYYWBGX2:4:11402:3004:17204\tchr3\t17376401\tchr4\t17467489\t+\t+"
          "\tUU\t60\t60");
    });
  };

  BENCHMARK_ADVANCED("sorting")
  (Catch::Benchmark::Chronometer meter) {
    const auto pixels = generate_pixels(1'000'000);
    std::vector<std::vector<Pixel<N>>> data(static_cast<std::size_t>(meter.runs()), pixels);
    meter.measure([&data](std::size_t i) {
      std::size_t num_ops{};
      std::sort(data[i].begin(), data[i].end(), [&](const auto& pixel1, const auto& pixel2) {
        ++num_ops;
        return pixel1 < pixel2;
      });
      return num_ops;
    });
  };
}

TEST_CASE("ThinPixel") {
  BENCHMARK_ADVANCED("from_coo w/table (uint32)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure([&bin_table]() {
      return Pixel<std::uint32_t>::from_coo(bin_table, "123456\t234567\t123");
    });
  };

  BENCHMARK_ADVANCED("from_coo w/table (double)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTable bin_table{hg38.begin(), hg38.end(), 1'000};
    meter.measure(
        [&bin_table]() { return Pixel<double>::from_coo(bin_table, "123456\t234567\t123.4567"); });
  };

  BENCHMARK_ADVANCED("from_coo wo/table (uint32)")
  (Catch::Benchmark::Chronometer meter) {
    meter.measure([]() { return ThinPixel<std::uint32_t>::from_coo("123456\t234567\t123"); });
  };

  BENCHMARK_ADVANCED("from_coo wo/table (double)")
  (Catch::Benchmark::Chronometer meter) {
    meter.measure([]() { return ThinPixel<double>::from_coo("123456\t234567\t123.4567"); });
  };

  BENCHMARK_ADVANCED("sorting")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table({hg38.begin(), hg38.end()}, 1'000);
    const auto pixels = generate_thin_pixels(1'000'000);
    std::vector<std::vector<ThinPixel<N>>> data(static_cast<std::size_t>(meter.runs()), pixels);
    meter.measure([&data](std::size_t i) {
      std::size_t num_ops{};
      std::sort(data[i].begin(), data[i].end(), [&](const auto& pixel1, const auto& pixel2) {
        ++num_ops;
        return pixel1 < pixel2;
      });
      return num_ops;
    });
  };
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)
