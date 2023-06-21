// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>
#include <readerwriterqueue/readerwriterqueue.h>
#include <spdlog/spdlog.h>

#include <future>

#include "hictk/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/tools/tools.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {

static bool missing_norm_or_interactions(const std::exception& e, hic::NormalizationMethod norm) {
  const std::string_view msg{e.what()};

  const auto missing_interactions =
      msg.find("unable to read file offset") != std::string_view::npos;

  const auto missing_norm_vect =
      msg.find(fmt::format(FMT_STRING("unable to find {} normalization vector"), norm)) !=
      std::string_view::npos;

  return missing_interactions || missing_norm_vect;
}

bool check_if_norm_exists(hic::HiCFile& f, hic::NormalizationMethod norm) {
  for (const auto& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    try {
      std::ignore = f.fetch(chrom.name(), norm);
      return true;
    } catch (const std::exception& e) {
      if (!missing_norm_or_interactions(e, norm)) {
        throw;
      }
    }
  }
  return false;
}

static std::vector<double> read_weights(hic::HiCFile& f, const BinTable& bins,
                                        hic::NormalizationMethod norm) {
  std::vector<double> weights{};
  weights.reserve(bins.size());
  std::size_t missing_norms = 0;
  for (const auto& chrom : bins.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto expected_length = (chrom.size() + bins.bin_size() - 1) / bins.bin_size();
    try {
      const auto weights_ = f.fetch(chrom.name(), norm).chrom1_norm();
      if (weights_.size() != expected_length) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("{} normalization vector for {} appears to be corrupted: "
                                   "expected {} values, found {}"),
                        norm, chrom.name(), expected_length, weights_.size()));
      }

      weights.insert(weights.end(), weights_.begin(), weights_.end());
      missing_norms = false;
    } catch (const std::exception& e) {
      if (!missing_norm_or_interactions(e, norm)) {
        throw;
      }
      weights.resize(weights.size() + expected_length, std::numeric_limits<double>::quiet_NaN());
      ++missing_norms;
    }
  }
  if (missing_norms == f.chromosomes().size() - 1) {
    spdlog::warn(FMT_STRING("[{}] {} normalization vector is missing. Filling "
                            "normalization vector with NaNs."),
                 bins.bin_size(), norm);
  }

  assert(weights.size() == bins.size());
  return weights;
}

template <typename CoolerFile>
static void copy_weights(hic::HiCFile& hf, CoolerFile& cf, hic::NormalizationMethod norm,
                         bool fail_if_missing) {
  if (norm == hic::NormalizationMethod::NONE) {
    return;
  }
  const auto dset_name = fmt::to_string(norm);

  const auto norm_exists = check_if_norm_exists(hf, norm);

  if (!norm_exists) {
    if (fail_if_missing) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find {} normalization vector for resolution {}"), norm,
                      hf.bins().bin_size()));
    }

    spdlog::warn(FMT_STRING("[{}] {} normalization vector is missing. SKIPPING!"),
                 hf.bins().bin_size(), norm);
    return;
  }

  spdlog::info(FMT_STRING("[{}] Processing {} normalization vector..."), hf.bins().bin_size(),
               norm);

  const auto weights = read_weights(hf, hf.bins(), norm);
  using T = std::remove_reference_t<decltype(cf)>;
  if constexpr (std::is_same_v<T, cooler::File>) {
    cf.write_weights(dset_name, weights.begin(), weights.end(), false, true);
  } else {
    cooler::File::write_weights(cf, dset_name, weights.begin(), weights.end(), false, true);
  }
}

[[nodiscard]] static cooler::File init_cooler(std::string_view uri, std::uint32_t resolution,
                                              std::string_view genome, const Reference& chroms) {
  auto attrs = cooler::StandardAttributes::init(resolution);
  attrs.assembly = genome.empty() ? "unknown" : std::string{genome};
  attrs.generated_by = fmt::format(FMT_STRING("hictk v{}"), hictk::config::version::str());

  return cooler::File::create_new_cooler(uri, chroms, resolution, true, attrs);
}

static Reference generate_reference(const std::filesystem::path& p, std::uint32_t res) {
  hic::HiCFile const hf(p, res);
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

template <typename N>
static void enqueue_pixels(const hic::HiCFile& hf,
                           moodycamel::BlockingReaderWriterQueue<Pixel<N>>& queue, bool quiet,
                           std::atomic<bool>& early_return,
                           std::size_t update_frequency = 10'000'000) {
  try {
    if (quiet) {
      update_frequency = (std::numeric_limits<std::size_t>::max)();
    }

    auto sel = hf.fetch();
    auto first = sel.begin<N>();
    auto last = sel.end<N>();

    auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; first != last && !early_return; ++i) {
      while (!queue.try_enqueue(*first)) {
        if (early_return) {
          return;
        }
        std::this_thread::sleep_for(std::chrono::microseconds(100));
      }
      ++first;

      if (i == update_frequency) {
        const auto t1 = std::chrono::steady_clock::now();
        const auto delta =
            static_cast<double>(
                std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
            1000.0;
        spdlog::info(
            FMT_STRING("[{}] Processing {:ucsc} at {:.0f} pixels/s (cache hit rate {:.2f}%)..."),
            hf.resolution(), first->coords.bin1, double(update_frequency) / delta,
            hf.block_cache_hit_rate() * 100);
        hf.reset_cache_stats();
        t0 = t1;
        i = 0;
      }
    }
    queue.enqueue(Pixel<N>{});
  } catch (...) {
    early_return = true;
    throw;
  }
}

template <typename N>
static std::size_t append_pixels(cooler::File& clr,
                                 moodycamel::BlockingReaderWriterQueue<Pixel<N>>& queue,
                                 std::atomic<bool>& early_return,
                                 std::size_t buffer_capacity = 100'000) {
  try {
    std::vector<Pixel<N>> buffer{buffer_capacity};
    buffer.clear();

    Pixel<N> value{};
    std::size_t nnz = 0;

    while (!early_return) {
      while (!queue.wait_dequeue_timed(value, std::chrono::milliseconds(25))) {
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

template <typename N>
static void convert_resolution_multi_threaded(
    hic::HiCFile& hf, std::string_view cooler_uri, const Reference& chromosomes,
    std::uint32_t resolution, std::string_view genome,
    const std::vector<hic::NormalizationMethod>& normalization_methods, bool fail_if_norm_not_found,
    bool quiet) {
  const auto t0 = std::chrono::steady_clock::now();

  spdlog::info(FMT_STRING("[{}] Begin processing {}bp matrix..."), hf.resolution(),
               hf.resolution());

  spdlog::debug(FMT_STRING("[{}] Block cache capacity: {}"), hf.resolution(), hf.cache_capacity());

  std::atomic<bool> early_return = false;
  moodycamel::BlockingReaderWriterQueue<Pixel<N>> queue{100'000};

  auto producer_fx = [&]() { return enqueue_pixels<N>(hf, queue, quiet, early_return); };
  auto producer = std::async(std::launch::async, producer_fx);

  auto clr = init_cooler(cooler_uri, resolution, genome, chromosomes);
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

  for (const auto norm : normalization_methods) {
    copy_weights(hf, clr, norm, fail_if_norm_not_found);
  }

  const auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  spdlog::info(FMT_STRING("[{}] DONE! Processed {} pixels across {} chromosomes in {:.2f}s"),
               resolution, nnz, hf.chromosomes().size() - 1, delta);
}

void convert_subcmd(const ConvertConfig& c) {
  assert(c.resolutions.size() > 0);

  assert(spdlog::default_logger());
  const auto t0 = std::chrono::steady_clock::now();

  if (c.resolutions.size() > 1) {
    cooler::init_mcool(c.output_cooler.string(), c.resolutions.begin(), c.resolutions.end(),
                       c.force);
  }

  const auto chroms = generate_reference(c.input_hic.string(), c.resolutions.front());
  hic::HiCFile hf(c.input_hic.string(), c.resolutions.front(), hic::MatrixType::observed,
                  hic::MatrixUnit::BP, c.block_cache_size);

  if (c.resolutions.size() == 1) {
    convert_resolution_multi_threaded<std::int32_t>(
        hf, c.output_cooler.string(), chroms, c.resolutions.front(), c.genome,
        c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable, c.quiet);
    return;
  }

  std::for_each(c.resolutions.rbegin(), c.resolutions.rend(), [&](const auto res) {
    hf.open(res);
    convert_resolution_multi_threaded<std::int32_t>(
        hf, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_cooler.string(), res), chroms,
        res, c.genome, c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable,
        c.quiet);
    hf.clear_cache();
  });

  const auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  spdlog::info(FMT_STRING("DONE! Processed {} resolution(s) in {:.2f}s!"), c.resolutions.size(),
               delta);
  spdlog::info(FMT_STRING("{} size: {:.2f} MB"), c.input_hic,
               static_cast<double>(std::filesystem::file_size(c.input_hic)) / 1.0e6);
  spdlog::info(FMT_STRING("{} size: {:.2f} MB"), c.output_cooler,
               static_cast<double>(std::filesystem::file_size(c.output_cooler)) / 1.0e6);
}
}  // namespace hictk::tools
