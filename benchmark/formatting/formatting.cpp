// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/benchmark/catch_constructor.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "./common.hpp"
#include "hictk/bin.hpp"
#include "hictk/bin_table_fixed.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/fmt.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/pixel.hpp"

using namespace hictk;

[[nodiscard]] static auto open_dev_null() {
  struct file_deleter {
    void operator()(std::FILE* fp) {
      std::fclose(fp);  // NOLINT
    }
  };

#ifdef _MSC_VER
  // NOLINTNEXTLINE(*-owning-memory)
  return std::unique_ptr<std::FILE, file_deleter>{std::fopen("nul", "w")};
#else
  // NOLINTNEXTLINE(*-owning-memory)
  return std::unique_ptr<std::FILE, file_deleter>{std::fopen("/dev/null", "w")};
#endif
}

template <typename N>
[[nodiscard]] static std::vector<Chromosome> to_chromosomes(
    const BinTableFixed& bin_table, const std::vector<ThinPixel<N>>& pixels) {
  std::vector<Chromosome> chroms{};
  chroms.reserve(pixels.size() * 2);

  for (const auto& p : pixels) {
    chroms.emplace_back(bin_table.at(p.bin1_id).chrom());
    chroms.emplace_back(bin_table.at(p.bin2_id).chrom());
  }

  return chroms;
}

template <typename N>
[[nodiscard]] static std::vector<Bin> to_bins(const BinTableFixed& bin_table,
                                              const std::vector<ThinPixel<N>>& pixels) {
  std::vector<Bin> bins{};
  bins.reserve(pixels.size() * 2);

  for (const auto& p : pixels) {
    bins.emplace_back(bin_table.at(p.bin1_id));
    bins.emplace_back(bin_table.at(p.bin2_id));
  }

  return bins;
}

template <typename N>
[[nodiscard]] static std::vector<GenomicInterval> to_genomic_intervals(
    const BinTableFixed& bin_table, const std::vector<ThinPixel<N>>& pixels) {
  std::vector<GenomicInterval> gis{};
  gis.reserve(pixels.size() * 2);

  for (const auto& p : pixels) {
    gis.emplace_back(bin_table.at(p.bin1_id).interval());
    gis.emplace_back(bin_table.at(p.bin2_id).interval());
  }

  return gis;
}

template <typename N>
[[nodiscard]] static std::vector<Pixel<N>> to_pixels(const BinTableFixed& bin_table,
                                                     const std::vector<ThinPixel<N>>& thin_pixels) {
  std::vector<Pixel<N>> pixels(thin_pixels.size());

  std::transform(thin_pixels.begin(), thin_pixels.end(), pixels.begin(),
                 [&](const ThinPixel<N>& tp) {
                   return Pixel<N>{bin_table.at(tp.bin1_id), bin_table.at(tp.bin2_id), tp.count};
                 });

  return pixels;
}

template <typename T>
[[nodiscard]] static std::vector<T> random_sample_with_replacement(const std::vector<T>& src,
                                                                   std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<T> dest{};
  while (dest.size() < size) {
    dest.insert(dest.end(), src.begin(), src.end());
  }

  dest.resize(size);
  std::shuffle(dest.begin(), dest.end(), rand_eng);

  return dest;
}

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Formatting Chromosome") {
  BENCHMARK_ADVANCED("wo/ compilation (TSV)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto chroms =
        random_sample_with_replacement(to_chromosomes(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&chroms, &fp]() {
      for (const auto& c : chroms) {
        fmt::print(fp.get(), FMT_STRING("{:tsv}"), c);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto chroms =
        random_sample_with_replacement(to_chromosomes(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&chroms, &fp]() {
      for (const auto& c : chroms) {
        fmt::print(fp.get(), FMT_STRING("{:ucsc}"), c);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (TSV)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto chroms =
        random_sample_with_replacement(to_chromosomes(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&chroms, &fp]() {
      for (const auto& c : chroms) {
        fmt::print(fp.get(), FMT_COMPILE("{:tsv}"), c);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto chroms =
        random_sample_with_replacement(to_chromosomes(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&chroms, &fp]() {
      for (const auto& c : chroms) {
        fmt::print(fp.get(), FMT_COMPILE("{:ucsc}"), c);
      }
    });
  };
}

TEST_CASE("Formatting GenomicInterval") {
  BENCHMARK_ADVANCED("wo/ compilation (BED)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto intervals =
        random_sample_with_replacement(to_genomic_intervals(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&intervals, &fp]() {
      for (const auto& gi : intervals) {
        fmt::print(fp.get(), FMT_STRING("{:bed}"), gi);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto intervals =
        random_sample_with_replacement(to_genomic_intervals(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&intervals, &fp]() {
      for (const auto& gi : intervals) {
        fmt::print(fp.get(), FMT_STRING("{:ucsc}"), gi);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (BED)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto intervals =
        random_sample_with_replacement(to_genomic_intervals(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&intervals, &fp]() {
      for (const auto& gi : intervals) {
        fmt::print(fp.get(), FMT_COMPILE("{:bed}"), gi);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto intervals =
        random_sample_with_replacement(to_genomic_intervals(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&intervals, &fp]() {
      for (const auto& gi : intervals) {
        fmt::print(fp.get(), FMT_COMPILE("{:ucsc}"), gi);
      }
    });
  };
}

TEST_CASE("Formatting Bin") {
  BENCHMARK_ADVANCED("wo/ compilation (raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_STRING("{:raw}"), bin);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (BED)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_STRING("{:bed}"), bin);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_STRING("{:ucsc}"), bin);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_COMPILE("{:raw}"), bin);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (BED)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_COMPILE("{:bed}"), bin);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (UCSC)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto bins = random_sample_with_replacement(to_bins(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&bins, &fp]() {
      for (const auto& bin : bins) {
        fmt::print(fp.get(), FMT_COMPILE("{:ucsc}"), bin);
      }
    });
  };
}

TEST_CASE("Formatting ThinPixel") {
  BENCHMARK_ADVANCED("wo/ compilation (int)")
  (Catch::Benchmark::Chronometer meter) {
    const auto pixels = random_sample_with_replacement(pixels_int, 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& tp : pixels) {
        fmt::print(fp.get(), FMT_STRING("{}"), tp);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (double)")
  (Catch::Benchmark::Chronometer meter) {
    const auto pixels = random_sample_with_replacement(pixels_fp, 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& tp : pixels) {
        fmt::print(fp.get(), FMT_STRING("{}"), tp);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (int)")
  (Catch::Benchmark::Chronometer meter) {
    const auto pixels = random_sample_with_replacement(pixels_int, 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& tp : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{}"), tp);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (double)")
  (Catch::Benchmark::Chronometer meter) {
    const auto pixels = random_sample_with_replacement(pixels_fp, 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& tp : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{}"), tp);
      }
    });
  };
}

TEST_CASE("Formatting Pixel") {
  BENCHMARK_ADVANCED("wo/ compilation (int; BG2)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_STRING("{:bg2}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (double; BG2)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_fp), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_STRING("{:bg2}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (int; BG2)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{:bg2}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (double; BG2)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_fp), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{:bg2}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (int; raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_STRING("{:raw}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("wo/ compilation (double; raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_fp), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_STRING("{:raw}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (int; raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_int), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{:raw}"), pxl);
      }
    });
  };

  BENCHMARK_ADVANCED("w/ compilation (double; raw)")
  (Catch::Benchmark::Chronometer meter) {
    const BinTableFixed bin_table{{hg38.begin(), hg38.end()}, 2'500'000};

    const auto pixels = random_sample_with_replacement(to_pixels(bin_table, pixels_fp), 100'000);
    auto fp = open_dev_null();

    meter.measure([&pixels, &fp]() {
      for (const auto& pxl : pixels) {
        fmt::print(fp.get(), FMT_COMPILE("{:raw}"), pxl);
      }
    });
  };
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)
