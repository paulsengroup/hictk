// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#if __has_include(<readerwriterqueue.h>)
#include <readerwriterqueue.h>
#else
#include <readerwriterqueue/readerwriterqueue.h>
#endif
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <future>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static bool missing_norm_or_interactions(const std::exception& e, balancing::Method norm) {
  const std::string_view msg{e.what()};

  const auto missing_interactions =
      msg.find("unable to read file offset") != std::string_view::npos;

  const auto missing_norm_vect =
      msg.find(fmt::format(FMT_STRING("unable to find {} normalization vector"), norm)) !=
      std::string_view::npos;

  return missing_interactions || missing_norm_vect;
}

static bool check_if_norm_exists(hic::File& f, balancing::Method norm) {
  return std::any_of(f.chromosomes().begin(), f.chromosomes().end(), [&](const Chromosome& chrom) {
    try {
      if (!chrom.is_all()) {
        std::ignore = f.fetch(chrom.name(), norm);
      }
      return true;
    } catch (const std::exception& e) {
      if (!missing_norm_or_interactions(e, norm)) {
        throw;
      }
    }
    return false;
  });
}

template <typename CoolerFile>
static void copy_weights(hic::File& hf, CoolerFile& cf, balancing::Method norm,
                         bool fail_if_missing) {
  if (norm == balancing::Method::NONE()) {
    return;
  }
  const auto dset_name = fmt::to_string(norm);

  const auto norm_exists = check_if_norm_exists(hf, norm);

  if (!norm_exists) {
    if (fail_if_missing) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find {} normalization vector for resolution {}"), norm,
                      hf.bins().resolution()));
    }

    SPDLOG_WARN(FMT_STRING("[{}] {} normalization vector is missing. SKIPPING!"),
                hf.bins().resolution(), norm);
    return;
  }

  SPDLOG_INFO(FMT_STRING("[{}] processing {} normalization vector..."), hf.bins().resolution(),
              norm);

  const auto weights = hf.normalization(norm)(balancing::Weights::Type::DIVISIVE);
  using T = std::remove_reference_t<decltype(cf)>;
  if constexpr (std::is_same_v<T, cooler::File>) {
    cf.write_weights(dset_name, weights.begin(), weights.end(), false, true);
  } else {
    cooler::File::write_weights(cf, dset_name, weights.begin(), weights.end(), false, true);
  }
}

template <typename PixelT>
[[nodiscard]] static cooler::File init_cooler(cooler::RootGroup entrypoint,
                                              std::uint32_t resolution, std::string_view genome,
                                              const Reference& chroms,
                                              std::uint32_t compression_lvl) {
  auto attrs = cooler::Attributes::init<PixelT>(resolution);
  attrs.assembly = genome.empty() ? "unknown" : std::string{genome};

  return cooler::File::create<PixelT>(std::move(entrypoint), chroms, resolution, attrs,
                                      cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);
}

template <typename PixelT>
[[nodiscard]] static cooler::File init_cooler(std::string_view uri, std::uint32_t resolution,
                                              std::string_view genome, const Reference& chroms,
                                              std::uint32_t compression_lvl) {
  auto attrs = cooler::Attributes::init<PixelT>(resolution);
  attrs.assembly = genome.empty() ? "unknown" : std::string{genome};

  return cooler::File::create<PixelT>(uri, chroms, resolution, true, attrs,
                                      cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);
}

static Reference generate_reference(const std::filesystem::path& p, std::uint32_t res) {
  hic::File const hf(p.string(), res);
  std::vector<std::string> names;
  std::vector<std::uint32_t> sizes;
  for (const auto& chrom : hf.chromosomes()) {
    if (!chrom.is_all()) {
      names.emplace_back(chrom.name());
      sizes.push_back(chrom.size());
    }
  }
  return {names.begin(), names.end(), sizes.begin()};
}

[[nodiscard]] static hic::PixelSelectorAll fetch_interactions_for_chromosome(
    hic::File& hf, const Chromosome& chrom1) {
  std::vector<hic::PixelSelector> selectors{};
  for (std::uint32_t chrom2_id = chrom1.id(); chrom2_id < hf.chromosomes().size(); ++chrom2_id) {
    const auto& chrom2 = hf.chromosomes().at(chrom2_id);
    if (chrom2.is_all()) {
      continue;
    }
    try {
      auto sel = hf.fetch(chrom1.name(), chrom2.name());
      if (!sel.empty()) {
        selectors.emplace_back(std::move(sel));
      }
    } catch (const std::exception& e) {
      const std::string_view msg{e.what()};
      const auto missing_norm = msg.find("unable to find") != std::string_view::npos &&
                                msg.find("normalization vector") != std::string_view::npos;
      if (!missing_norm) {
        throw;
      }
    }
  }

  return hic::PixelSelectorAll{std::move(selectors)};
}

template <typename N>
static void enqueue_pixels(hic::File& hf,
                           moodycamel::BlockingReaderWriterQueue<ThinPixel<N>>& queue,
                           std::atomic<bool>& early_return,
                           std::size_t update_frequency = 10'000'000) {
  try {
    std::size_t i = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (std::uint32_t chrom1_id = 0; chrom1_id < hf.chromosomes().size(); ++chrom1_id) {
      hf.purge_footer_cache();
      hf.clear_cache();

      const auto& chrom1 = hf.chromosomes().at(chrom1_id);
      if (chrom1.is_all()) {
        continue;
      }

      const auto sel = fetch_interactions_for_chromosome(hf, chrom1);
      auto first = sel.begin<N>();
      auto last = sel.end<N>();

      for (; first != last && !early_return; ++i) {
        while (!queue.try_enqueue(*first)) {
          if (early_return) {
            return;
          }
        }
        ++first;

        if (i == update_frequency) {
          const auto t1 = std::chrono::steady_clock::now();
          const auto delta =
              static_cast<double>(
                  std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
              1000.0;
          const auto bin1 = hf.bins().at(first->bin1_id);
          SPDLOG_INFO(
              FMT_STRING("[{}] processing {:ucsc} at {:.0f} pixels/s (cache hit rate {:.2f}%)..."),
              hf.resolution(), bin1, double(update_frequency) / delta,
              hf.block_cache_hit_rate() * 100);
          hf.reset_cache_stats();
          t0 = t1;
          i = 0;
        }
      }
    }
    queue.enqueue(ThinPixel<N>{});
  } catch (...) {
    early_return = true;
    throw;
  }
}

template <typename N>
static std::size_t append_pixels(cooler::File& clr,
                                 moodycamel::BlockingReaderWriterQueue<ThinPixel<N>>& queue,
                                 std::atomic<bool>& early_return,
                                 std::size_t buffer_capacity = 100'000) {
  try {
    std::vector<ThinPixel<N>> buffer{buffer_capacity};
    buffer.clear();

    ThinPixel<N> value{};
    std::size_t nnz = 0;

    while (!early_return) {
      while (!queue.try_dequeue(value)) {
        if (early_return) {
          return nnz;
        }
      }

      if (!value) {
        break;
      }

      buffer.push_back(value);
      ++nnz;

      if (buffer.size() == buffer.capacity()) {
        clr.append_pixels(buffer.begin(), buffer.end());
        buffer.clear();
      }
    }

    if (!buffer.empty()) {
      clr.append_pixels(buffer.begin(), buffer.end());
    }
    return nnz + buffer.size();
  } catch (...) {
    early_return = true;
    throw;
  }
}

template <typename N>  // NOLINTNEXTLINE(*-rvalue-reference-param-not-moved)
static void convert_resolution_multi_threaded(hic::File& hf, cooler::File&& clr,
                                              std::vector<balancing::Method> normalization_methods,
                                              bool fail_if_norm_not_found) {
  const auto t0 = std::chrono::steady_clock::now();

  if (normalization_methods.empty()) {
    normalization_methods = hf.avail_normalizations();
  }

  SPDLOG_INFO(FMT_STRING("[{}] begin processing {}bp matrix..."), hf.resolution(), hf.resolution());

  std::atomic<bool> early_return = false;
  moodycamel::BlockingReaderWriterQueue<ThinPixel<N>> queue{1'000'000};

  auto producer_fx = [&]() { return enqueue_pixels<N>(hf, queue, early_return); };
  auto producer = std::async(std::launch::async, producer_fx);

  auto consumer_fx = [&]() { return append_pixels<N>(clr, queue, early_return); };
  auto consumer = std::async(std::launch::async, consumer_fx);

  try {
    producer.get();
  } catch (const std::exception& e) {
    early_return = true;
    throw std::runtime_error(fmt::format(
        FMT_STRING("exception raised while reading interactions from input file: {}"), e.what()));
  }

  std::size_t nnz = 0;
  try {
    nnz = consumer.get();
  } catch (const std::exception& e) {
    early_return = true;
    throw std::runtime_error(fmt::format(
        FMT_STRING("exception raised while writing interactions to output file: {}"), e.what()));
  }

  for (const auto& norm : normalization_methods) {
    copy_weights(hf, clr, norm, fail_if_norm_not_found);
  }

  const auto resolution = clr.resolution();
  clr.close();
  const auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  SPDLOG_INFO(FMT_STRING("[{}] DONE! Processed {} pixels across {} chromosomes in {:.2f}s"),
              resolution, nnz, hf.chromosomes().size() - 1, delta);
}

[[nodiscard]] static std::variant<std::int32_t, float> infer_count_type(
    const std::filesystem::path& p, std::size_t max_sample_size = 1'000'000) {
  const auto base_resolution = hic::utils::list_resolutions(p, true).front();
  const hic::File f(p, base_resolution);
  const auto sel = f.fetch();

  auto first_pixel = sel.begin<float>();
  const auto last_pixel = sel.end<float>();

  for (std::size_t i = 0; i < max_sample_size && first_pixel != last_pixel; ++i, ++first_pixel) {
    const auto pixel = *first_pixel;
    if (pixel.count != std::floor(pixel.count)) {
      return {float{}};
    }
  }

  return {std::int32_t{}};
}

void hic_to_cool(const ConvertConfig& c) {
  assert(!c.resolutions.empty());

  std::variant<std::int32_t, float> count_type{std::int32_t{}};

  if (c.count_type == "auto") {
    count_type = infer_count_type(c.path_to_input);
  } else if (c.count_type == "int") {
    count_type = std::int32_t{};
  } else {
    assert(c.count_type == "float");
    count_type = float{};
  }

  std::visit(
      [&]([[maybe_unused]] auto count_type_) {
        using PixelT = decltype(count_type_);

        const auto chroms = generate_reference(c.path_to_input.string(), c.resolutions.front());
        hic::File hf(c.path_to_input.string(), c.resolutions.front());
        assert(spdlog::default_logger());

        if (c.resolutions.size() == 1) {
          convert_resolution_multi_threaded<PixelT>(
              hf,
              init_cooler<PixelT>(c.path_to_output.string(), c.resolutions.front(), c.genome,
                                  chroms, c.compression_lvl),
              c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable);
          return;
        }

        auto mclr = cooler::MultiResFile::create(c.path_to_output.string(), chroms, c.force);

        std::for_each(c.resolutions.begin(), c.resolutions.end(), [&](const auto res) {
          hf.open(res);
          auto attrs = cooler::Attributes::init<PixelT>(res);
          attrs.assembly = c.genome.empty() ? "unknown" : std::string{c.genome};
          convert_resolution_multi_threaded<PixelT>(
              hf,
              init_cooler<PixelT>(mclr.init_resolution(res), res, c.genome, chroms,
                                  c.compression_lvl),
              c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable);
          hf.clear_cache();
        });
      },
      count_type);
}
}  // namespace hictk::tools
