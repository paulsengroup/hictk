

// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/benchmark/utils.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/reference.hpp"

namespace hictk::benchmark {

template <typename N>
[[nodiscard]] constexpr std::string_view type_name() {
  static_assert(std::is_arithmetic_v<N>);
  static_assert(sizeof(float) == 4);   // NOLINT(*-avoid-magic-numbers)
  static_assert(sizeof(double) == 8);  // NOLINT(*-avoid-magic-numbers)
  if constexpr (std::is_same_v<N, std::uint8_t>) {
    return "std::uint8_t";
  }
  if constexpr (std::is_same_v<N, std::uint16_t>) {
    return "std::uint16_t";
  }
  if constexpr (std::is_same_v<N, std::uint32_t>) {
    return "std::uint32_t";
  }
  if constexpr (std::is_same_v<N, std::uint64_t>) {
    return "std::uint64_t";
  }

  if constexpr (std::is_same_v<N, std::int8_t>) {
    return "std::int8_t";
  }
  if constexpr (std::is_same_v<N, std::int16_t>) {
    return "std::int16_t";
  }
  if constexpr (std::is_same_v<N, std::int32_t>) {
    return "std::int32_t";
  }
  if constexpr (std::is_same_v<N, std::int64_t>) {
    return "std::int64_t";
  }

  if constexpr (std::is_same_v<N, float>) {
    return "std::float32_t";
  }
  if constexpr (std::is_same_v<N, double>) {
    return "std::float64_t";
  }

  throw std::logic_error("unsupported type");
}

struct Params {
  std::string name;
  std::string tags;
  std::filesystem::path path;
  std::uint32_t resolution{};
  std::string_view range1;
  std::string_view range2;
  balancing::Method normalization;
};

// This class is necessary to work around some limitations of Catch2.
// Basically what we are after are parametrized benchmarks with automatic name and tag generation.
// This is not possible out of the box as of Catch2 v3.7.1 for the following reasons:
// - We can only assign tags to test cases: that means we need to have one benchmark per test case.
// - Catch supports dynamic test registration through the REGISTER_TEST_CASE macro.
//   Unfortunately, this macro only supports registering functions with signature void foo();
//   which prevents us from passing benchmark params directly to the benchmark function.
// To work around these limitations I have implemented the TestCaseGenerator class (see below).
// This class can be constructed given one or more std::array with the params. This allows us to
// compute the number of parameter combinations at compile time. This is useful because we can use
// std::index_sequence in combination with expression folding to call template <std::size_t I> void
// run_benchmark(). This function then uses the template parameter to fetch the appropriate set of
// parameters from a TestCaseGenerator class that has been declared in the global namespace as a
// static const variable.
template <std::size_t S1, std::size_t S2 = 0, std::size_t S3 = 0, std::size_t S4 = 0,
          std::size_t S5 = 0>
class TestCaseGenerator {
 public:
  [[nodiscard]] static constexpr std::size_t size() noexcept {
    return S1 * std::max(std::size_t{1}, S2) * std::max(std::size_t{1}, S3) *
           std::max(std::size_t{1}, S4) * std::max(std::size_t{1}, S5);
  }

 private:
  std::vector<Params> _params{};
  static constexpr std::size_t _chunk_size{32};

  [[nodiscard]] static std::uint64_t compute_num_pixels_ub(std::string_view range1,
                                                           std::string_view range2,
                                                           std::uint32_t resolution) {
    const auto gi1 = GenomicInterval::parse_ucsc(std::string{range1});
    const auto gi2 = GenomicInterval::parse_ucsc(std::string{range2});

    const auto size1 = std::get<2>(gi1) - std::get<1>(gi1);
    const auto size2 = std::get<2>(gi2) - std::get<1>(gi2);

    const auto nbins1 = static_cast<std::uint64_t>((size1 + resolution - 1) / resolution);
    const auto nbins2 = static_cast<std::uint64_t>((size2 + resolution - 1) / resolution);

    return nbins1 * nbins2;
  }

  template <typename N>
  [[nodiscard]] static std::string generate_tags(const std::filesystem::path& path,
                                                 std::string_view range1, std::string_view range2,
                                                 std::uint32_t resolution) {
    const auto ext = path.extension();
    assert(!ext.empty());
    auto tags = fmt::format(FMT_STRING("[benchmark][interaction_fetching][{}][{}bp]"),
                            ext.string().substr(1), resolution);
    if (range1 == "GW") {
      assert(range2 == "GW");
      tags += "[gw]";
    } else if (range1 == range2) {
      tags += "[cis]";
    } else {
      tags += "[trans]";
    }

    if (range1 == "GW") {
      assert(range1 == "GW");
      tags += "[large]";
    } else {
      const auto num_pixels = compute_num_pixels_ub(range1, range2, resolution);
      if (num_pixels < 100'000) {  // NOLINT(*-avoid-magic-numbers)
        tags += "[small]";
      } else if (num_pixels < 2'500'000) {  // NOLINT(*-avoid-magic-numbers)
        tags += "[medium]";
      } else {
        tags += "[large]";
      }
    }

    tags += fmt::format(FMT_STRING("[{}]"), type_name<N>());

    return tags;
  }

 public:
  TestCaseGenerator() = delete;
  // NOLINTNEXTLINE(*-function-cognitive-complexity)
  TestCaseGenerator(std::string_view title, std::array<std::string_view, S1> files,
                    std::array<std::uint32_t, S2> resolutions,
                    std::array<std::string_view, S3> ranges1,
                    std::array<std::string_view, S4> ranges2,
                    std::array<balancing::Method, S5> normalizations)
      : _params(size()) {
    if constexpr (size() == 0) {
      throw std::logic_error("size cannot be 0");
    }
    std::size_t i = 0;
    const auto test_name = generate_test_name(title, false);
    for (const auto& f : files) {
      for (const auto& res : resolutions) {
        for (const auto& r1 : ranges1) {
          for (const auto& r2 : ranges2) {
            for (const auto& norm : normalizations) {
              const auto int_counts = norm == balancing::Method::NONE();
              std::filesystem::path path{f};

              const auto format = path.extension().string().substr(1);

              _params[i].name = fmt::format(FMT_STRING("{{"
                                                       "{}, "
                                                       "\"format\": \"{}\", "
                                                       "\"range1\": \"{}\", "
                                                       "\"range2\": \"{}\", "
                                                       "\"resolution\": {}, "
                                                       "\"sorted\": true, "
                                                       "\"count-type\": \"{}\""
                                                       "}}"),
                                            test_name, format, r1, r2, res,
                                            int_counts ? "std::uint32_t" : "std::float64_t", norm);
              _params[i].tags = int_counts ? generate_tags<std::uint32_t>(path, r1, r2, res)
                                           : generate_tags<double>(path, r1, r2, res);
              _params[i].path = std::move(path);
              _params[i].resolution = res;
              _params[i].range1 = r1;
              _params[i].range2 = r2;
              _params[i++].normalization = norm;
            }
          }
        }
      }
    }
  }

  [[nodiscard]] constexpr auto operator[](std::size_t i) const noexcept -> const Params& {
    assert(i < size());
    return _params[i];
  }

  [[nodiscard]] static constexpr std::size_t chunk_size() noexcept { return _chunk_size; }

  [[nodiscard]] static constexpr std::size_t num_chunks() noexcept {
    return (size() + chunk_size() - 1) / chunk_size();
  }
};

// This macro defines the boilerplate required to register one test case (which corresponds to one
// benchmark) for each parameter combination generated by an instance of the TestCaseGenerator
// defined above.
// In principle, we should be able to use fold expressions plus a std::index_sequence to call
// template <std::size_t I> void run_benchmark() to run a single test case.
// Unfortunately there are scenarios where a single TestCaseGenerator yields too many parameter
// combinations to be handled by a single fold expression.
// The boilerplate defined below basically splits a std::index_sequence into multiple chunks.
// Then submits one chunk at a time.
#define HICTK_REGISTER_BENCHMARKS(_ParamGenerator, _BenchmarkRunner)                  \
  namespace internal {                                                                \
  template <std::size_t I, std::size_t J>                                             \
  static void register_benchmark() {                                                  \
    constexpr auto IDX = (I * _ParamGenerator.chunk_size()) + J;                      \
    if constexpr (IDX < _ParamGenerator.size()) {                                     \
      REGISTER_TEST_CASE(_BenchmarkRunner<IDX>, _ParamGenerator[IDX].name,            \
                         _ParamGenerator[IDX].tags);                                  \
    }                                                                                 \
  }                                                                                   \
                                                                                      \
  template <std::size_t Chunk, std::size_t... Is>                                     \
  static void register_benchmarks_chunk(std::index_sequence<Is...>) {                 \
    (register_benchmark<Chunk, Is>(), ...);                                           \
  }                                                                                   \
                                                                                      \
  template <std::size_t... Chunks>                                                    \
  static void register_benchmarks(std::index_sequence<Chunks...>) {                   \
    constexpr auto CHUNK_SIZE = _ParamGenerator.chunk_size();                         \
    (register_benchmarks_chunk<Chunks>(std::make_index_sequence<CHUNK_SIZE>{}), ...); \
  }                                                                                   \
  }                                                                                   \
                                                                                      \
  static void register_benchmarks() {                                                 \
    constexpr auto NUM_CHUNKS = _ParamGenerator.num_chunks();                         \
    internal::register_benchmarks(std::make_index_sequence<NUM_CHUNKS>{});            \
  }

[[nodiscard]] inline std::pair<std::string, std::string> generate_query(
    std::mt19937_64& rand_eng, const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
    double avg_height, double avg_width, double height_std, double width_std) {
  assert(chrom1 <= chrom2);

  const auto pos1 = std::uniform_int_distribution<std::uint32_t>{0U, chrom1.size() - 1}(rand_eng);
  const auto pos2 = std::uniform_int_distribution<std::uint32_t>{0U, chrom2.size() - 1}(rand_eng);

  const auto height = static_cast<std::uint32_t>(
      std::clamp(std::normal_distribution<double>{avg_height, height_std}(rand_eng), 1.0,
                 static_cast<double>(chrom1.size())));
  const auto width = static_cast<std::uint32_t>(
      std::clamp(std::normal_distribution<double>{avg_width, width_std}(rand_eng), 1.0,
                 static_cast<double>(chrom2.size())));

  auto start1 = height >= pos1 ? std::uint32_t{0} : pos1 - height;
  auto start2 = width >= pos2 ? std::uint32_t{0} : pos2 - width;

  if (chrom1 == chrom2 && start1 > start2) {
    std::swap(start1, start2);
  }

  auto end1 = std::min(start1 + height, chrom1.size());
  auto end2 = std::min(start2 + width, chrom2.size());

  return std::make_pair(fmt::format(FMT_STRING("{}:{}-{}"), chrom1.name(), start1, end1),
                        fmt::format(FMT_STRING("{}:{}-{}"), chrom2.name(), start2, end2));
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries(
    const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2, std::size_t num_queries,
    double avg_height, double avg_width, double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    return generate_query(rand_eng, chrom1, chrom2, avg_height, avg_width, height_std, width_std);
  });

  return queries;
}

[[nodiscard]] inline std::discrete_distribution<std::uint32_t> init_chromosome_selector(
    const hictk::Reference& chroms) {
  std::vector<std::uint32_t> weights{};
  for (const auto& chrom : chroms) {
    if (chrom.is_all()) {
      weights.push_back(0);
      continue;
    }
    weights.push_back(chrom.size());
  }

  return {weights.begin(), weights.end()};
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries_cis(
    const hictk::Reference& chroms, std::size_t num_queries, double avg_height, double avg_width,
    double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  auto chrom_selector = init_chromosome_selector(chroms);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    const auto& chrom = chroms.at(chrom_selector(rand_eng));
    return generate_query(rand_eng, chrom, chrom, avg_height, avg_width, height_std, width_std);
  });

  return queries;
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries_trans(
    const hictk::Reference& chroms, std::size_t num_queries, double avg_height, double avg_width,
    double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  auto chrom_selector = init_chromosome_selector(chroms);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    const auto& chrom1 = chroms.at(chrom_selector(rand_eng));
    while (true) {
      const auto& chrom2 = chroms.at(chrom_selector(rand_eng));
      if (chrom1 == chrom2) {
        continue;
      }
      return generate_query(rand_eng, chrom1, chrom2, avg_height, avg_width, height_std, width_std);
    }
  });

  return queries;
}

template <typename N, typename File>
[[nodiscard]] inline std::ptrdiff_t count_nnz(const File& file, std::string_view range1,
                                              std::string_view range2,
                                              const hictk::balancing::Method& normalization) {
  const auto sel = file.fetch(range1, range2, normalization);

  return std::distance(sel.template begin<N>(), sel.template end<N>());
}

template <typename N, typename File>
[[nodiscard]] inline std::ptrdiff_t count_nnz(const File& file, std::size_t max_num_pixels,
                                              const hictk::balancing::Method& normalization) {
  const auto sel = file.fetch(normalization);
  auto first = sel.template begin<N>();
  auto last = sel.template end<N>();

  std::ptrdiff_t i{};
  // clang-format off
  while (++first != last && ++i != static_cast<std::ptrdiff_t>(max_num_pixels));  // NOLINT
  // clang-format on
  return i;
}

template <typename N, typename File>
[[nodiscard]] inline std::ptrdiff_t count_nnz_unsorted(
    const File& file, std::string_view range1, std::string_view range2,
    const hictk::balancing::Method& normalization) {
  const auto sel = file.fetch(range1, range2, normalization);

  return std::distance(sel.template begin<N>(false), sel.template end<N>());
}

template <typename N, typename File>
[[nodiscard]] inline std::ptrdiff_t count_nnz_unsorted(
    const File& file, std::size_t max_num_pixels, const hictk::balancing::Method& normalization) {
  const auto sel = file.fetch(normalization);
  auto first = sel.template begin<N>(false);
  auto last = sel.template end<N>();

  std::ptrdiff_t i{};
  // clang-format off
  while (++first != last && ++i != static_cast<std::ptrdiff_t>(max_num_pixels));  // NOLINT
  // clang-format on
  return i;
}

}  // namespace hictk::benchmark
