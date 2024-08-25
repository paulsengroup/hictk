// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <arrow/table.h>
#include <fmt/format.h>
#include <pybind11/embed.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

#include "hictk/file.hpp"
#include "hictk/fuzzer/config.hpp"
#include "hictk/fuzzer/cooler.hpp"
#include "hictk/fuzzer/tools.hpp"
#include "hictk/fuzzer/validators.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/to_dataframe.hpp"

namespace hictk::fuzzer {

// clang-format off
using PixelBuffer =
  std::variant<
      std::vector<ThinPixel<std::int32_t>>,
      std::vector<ThinPixel<double>>,
      std::vector<Pixel<std::int32_t>>,
      std::vector<Pixel<double>>>;
// clang-format on

template <typename N>
static void to_vector(std::vector<ThinPixel<N>>& buff, const std::shared_ptr<arrow::Table>& data) {
  buff.clear();
  buff.reserve(static_cast<std::size_t>(data->num_rows()));

  const auto bin1_id_chunks = data->GetColumnByName("bin1_id")->chunks();
  const auto bin2_id_chunks = data->GetColumnByName("bin2_id")->chunks();
  const auto count_chunks = data->GetColumnByName("count")->chunks();

  for (std::size_t chunk_id = 0; chunk_id < bin1_id_chunks.size(); ++chunk_id) {
    const auto& bin1_id_buffers = bin1_id_chunks[chunk_id]->data()->buffers;
    const auto& bin2_id_buffers = bin2_id_chunks[chunk_id]->data()->buffers;
    const auto& count_buffers = count_chunks[chunk_id]->data()->buffers;

    for (std::size_t buffer_id = 0; buffer_id < bin1_id_buffers.size(); ++buffer_id) {
      if (!bin1_id_buffers[buffer_id]) {
        continue;
      }
      const auto size =
          bin1_id_buffers[buffer_id]->size() / static_cast<std::int64_t>(sizeof(std::uint64_t));
      const auto* bin1_ids = bin1_id_buffers[buffer_id]->data_as<std::uint64_t>();
      const auto* bin2_ids = bin2_id_buffers[buffer_id]->data_as<std::uint64_t>();
      const auto* counts = count_buffers[buffer_id]->data_as<N>();
      assert(bin1_ids);
      assert(bin2_ids);
      assert(counts);
      for (std::int64_t i = 0; i < size; ++i) {
        buff.emplace_back(ThinPixel<N>{*(bin1_ids + i), *(bin2_ids + i), *(counts + i)});
      }
    }
  }
}

template <typename N>
static void to_vector(const Reference& chroms, std::vector<Pixel<N>>& buff,
                      const std::shared_ptr<arrow::Table>& data) {
  buff.clear();
  buff.reserve(static_cast<std::size_t>(data->num_rows()));

  const auto chrom1_id_chunks = data->GetColumnByName("chrom1")->chunks();
  const auto start1_chunks = data->GetColumnByName("start1")->chunks();
  const auto end1_chunks = data->GetColumnByName("end1")->chunks();
  const auto chrom2_id_chunks = data->GetColumnByName("chrom2")->chunks();
  const auto start2_chunks = data->GetColumnByName("start2")->chunks();
  const auto end2_chunks = data->GetColumnByName("end2")->chunks();
  const auto count_chunks = data->GetColumnByName("count")->chunks();

  for (std::size_t chunk_id = 0; chunk_id < chrom1_id_chunks.size(); ++chunk_id) {
    const auto& chrom1_id_buffers = chrom1_id_chunks[chunk_id]->data()->buffers;
    const auto& start1_buffers = start1_chunks[chunk_id]->data()->buffers;
    const auto& end1_buffers = end1_chunks[chunk_id]->data()->buffers;
    const auto& chrom2_id_buffers = chrom2_id_chunks[chunk_id]->data()->buffers;
    const auto& start2_buffers = start2_chunks[chunk_id]->data()->buffers;
    const auto& end2_buffers = end2_chunks[chunk_id]->data()->buffers;
    const auto& count_buffers = count_chunks[chunk_id]->data()->buffers;

    for (std::size_t buffer_id = 0; buffer_id < chrom1_id_buffers.size(); ++buffer_id) {
      if (!chrom1_id_buffers[buffer_id]) {
        continue;
      }
      const auto size =
          chrom1_id_buffers[buffer_id]->size() / static_cast<std::int64_t>(sizeof(std::uint32_t));
      const auto* chrom1_id = chrom1_id_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* start1 = start1_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* end1 = end1_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* chrom2_id = chrom2_id_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* start2 = start2_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* end2 = end2_buffers[buffer_id]->data_as<std::uint32_t>();
      const auto* counts = count_buffers[buffer_id]->data_as<N>();
      assert(chrom1_id);
      assert(start1);
      assert(end1);
      assert(chrom2_id);
      assert(start2);
      assert(end2);
      assert(counts);
      for (std::int64_t i = 0; i < size; ++i) {
        buff.emplace_back(Pixel{chroms.at(*(chrom1_id + i)), *(start1 + i), *(end1 + i),
                                chroms.at(*(chrom2_id + i)), *(start2 + i), *(end2 + i),
                                *(counts + i)});
      }
    }
  }
}

struct Query {
  Chromosome chrom{};
  double start_pos{};
  double end_pos{};

  [[nodiscard]] std::string to_string() const {
    return fmt::format(FMT_STRING("{}:{:.0f}-{:.0f}"), chrom.name(), start_pos, end_pos);
  }
};

[[nodiscard]] static Query generate_query_1d(
    const hictk::Reference& chroms, std::mt19937_64& rand_eng,
    std::discrete_distribution<std::uint32_t>& chrom_sampler, double mean_length,
    double stddev_length) {
  assert(!chroms.empty());

  const auto& chrom = chroms[chrom_sampler(rand_eng)];
  const auto query_length = static_cast<std::uint32_t>(
      std::normal_distribution<double>{mean_length, stddev_length}(rand_eng));

  const auto center_pos =
      std::uniform_real_distribution<double>{0.0, static_cast<double>(chrom.size())}(rand_eng);
  const auto start_pos = std::max(0.0, center_pos - (query_length / 2));
  const auto end_pos = std::min(static_cast<double>(chrom.size()), start_pos + query_length);

  return {chrom, start_pos, end_pos};
}

[[nodiscard]] static std::pair<Query, Query> generate_query_2d(
    const hictk::Reference& chroms, std::mt19937_64& rand_eng,
    std::discrete_distribution<std::uint32_t>& chrom_sampler, double mean_length,
    double stddev_length) {
  auto q1 = generate_query_1d(chroms, rand_eng, chrom_sampler, mean_length, stddev_length);
  auto q2 = generate_query_1d(chroms, rand_eng, chrom_sampler, mean_length, stddev_length);

  if (q1.chrom.id() > q2.chrom.id()) {
    std::swap(q1, q2);
  }

  if (q1.chrom == q2.chrom && q1.start_pos > q2.start_pos) {
    std::swap(q1, q2);
  }

  return std::make_pair(q1, q2);
}

[[nodiscard]] std::discrete_distribution<std::uint32_t> init_chrom_sampler(
    const hictk::Reference& chroms) {
  std::uint64_t genome_size{};
  std::vector<double> weights(chroms.size());
  std::transform(chroms.begin(), chroms.end(), weights.begin(),
                 [&](const hictk::Chromosome& chrom) {
                   genome_size += chrom.size();
                   return static_cast<double>(chrom.size());
                 });
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [&](const auto n) { return n / static_cast<double>(genome_size); });

  return {weights.begin(), weights.end()};
}

[[nodiscard]] static std::chrono::seconds compute_elapsed_time(
    const std::chrono::system_clock::time_point& t0) {
  const auto delta = std::chrono::system_clock::now() - t0;
  return std::chrono::duration_cast<std::chrono::seconds>(delta);
}

template <typename PixelT>
static void fetch(const Reference& chroms, cooler::Cooler& clr, std::string_view range1,
                  std::string_view range2, std::string_view normalization,
                  std::vector<PixelT>& buffer) {
  using N = decltype(std::declval<PixelT>().count);

  if constexpr (std::is_same_v<remove_cvref_t<PixelT>, Pixel<N>>) {
    cooler::BG2DataFrame<N> df{};
    clr.fetch_df(df, range1, range2, normalization);
    df.to_vector(chroms, buffer);
  } else {
    cooler::COODataFrame<N> df{};
    clr.fetch_df(df, range1, range2, normalization);
    df.to_vector(buffer);
  }
  assert(std::is_sorted(buffer.begin(), buffer.end()));
}

template <typename PixelT>
static void fetch(const hictk::File& f, std::string_view range1, std::string_view range2,
                  std::string_view normalization, std::vector<PixelT>& buffer) {
  using N = decltype(std::declval<PixelT>().count);

  const auto sel = f.fetch(range1, range2, balancing::Method{normalization});

  if constexpr (std::is_same_v<remove_cvref_t<PixelT>, ThinPixel<N>>) {
    to_vector(buffer, transformers::ToDataFrame(sel.begin<N>(), sel.end<N>(),
                                                transformers::DataFrameFormat::COO)());
  } else {
    to_vector(f.chromosomes(), buffer,
              transformers::ToDataFrame(sel.begin<N>(), sel.end<N>(),
                                        transformers::DataFrameFormat::BG2, f.bins_ptr())());
  }
}

[[nodiscard]] PixelBuffer init_pixel_buffer(const Config& c) {
  const auto int_count = c.normalization.empty() || c.normalization == "NONE";
  const auto thin_pixel = !c.join;

  if (thin_pixel && int_count) {
    return std::vector<ThinPixel<std::int32_t>>{};
  }
  if (thin_pixel && !int_count) {
    return std::vector<ThinPixel<double>>{};
  }

  if (int_count) {
    return std::vector<Pixel<std::int32_t>>{};
  }
  return std::vector<Pixel<double>>{};
}

int launch_worker_subcommand(const Config& c) {
  [[maybe_unused]] const pybind11::scoped_interpreter guard{};

  std::size_t num_tests{};
  std::size_t num_failures{};
  const std::chrono::microseconds duration{static_cast<std::int64_t>(c.duration * 1.0e6)};

  try {
    spdlog::info(FMT_STRING("[task_id] seed: {}"), c.seed);
    std::mt19937_64 rand_eng{c.seed};

    const hictk::File tgt(c.reference_uri, c.resolution);
    cooler::Cooler ref(c.resolution == 0 ? c.reference_uri.string()
                                         : fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                                       c.reference_uri.string(), c.resolution));

    const auto chroms = tgt.chromosomes().remove_ALL();
    auto chrom_sampler = init_chrom_sampler(chroms);

    auto expected_buffer = init_pixel_buffer(c);
    auto found_buffer = init_pixel_buffer(c);

    const auto t0 = std::chrono::system_clock::now();

    return std::visit(
        [&](auto& expected) {
          using BufferT = remove_cvref_t<decltype(expected)>;
          auto& found = std::get<BufferT>(found_buffer);

          while (compute_elapsed_time(t0) < duration) {
            const auto [q1, q2] = generate_query_2d(chroms, rand_eng, chrom_sampler,
                                                    c.query_length_avg, c.query_length_std);
            const auto range1 = q1.to_string();
            const auto range2 = q2.to_string();

            fetch(tgt.chromosomes(), ref, range1, range2, c.normalization, expected);
            fetch(tgt, range1, range2, c.normalization, found);

            ++num_tests;
            num_failures += !compare_pixels(range1, range2, expected, found);  // NOLINT
          }

          const auto num_successes = num_tests - num_failures;
          const auto failure_ratio =
              100.0 * static_cast<double>(num_successes) / static_cast<double>(num_tests);

          if (num_failures == 0) {
            spdlog::info(FMT_STRING("[task_id] Score: {:.4g} ({} success and {} failures)."),
                         failure_ratio, num_successes, num_failures);
            return 0;
          }
          spdlog::warn(FMT_STRING("[task_id] Score: {:.4g} ({} success and {} failures)."),
                       failure_ratio, num_successes, num_failures);
          return 1;
        },
        expected_buffer);
  } catch (const std::exception& e) {
    // wrap python exceptions while we still hold the scoped_interpreter guard
    throw std::runtime_error(e.what());
  }
}

}  // namespace hictk::fuzzer
