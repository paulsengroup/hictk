// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

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
    } catch (const std::exception& e) {
      if (!missing_norm_or_interactions(e, norm)) {
        throw;
      }
      spdlog::warn(FMT_STRING("[{}] {} normalization vector for {} is missing. Filling "
                              "normalization vector with NaNs."),
                   bins.bin_size(), norm, chrom.name());
      weights.resize(weights.size() + expected_length, std::numeric_limits<double>::quiet_NaN());
    }
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

static void convert_resolution(hic::HiCFile& hf, std::string_view cooler_uri,
                               const Reference& chromosomes, std::uint32_t resolution,
                               std::string_view genome,
                               const std::vector<hic::NormalizationMethod>& normalization_methods,
                               bool fail_if_norm_not_found) {
  const auto t0 = std::chrono::steady_clock::now();
  std::size_t nnz = 0;
  auto cf = init_cooler(cooler_uri, resolution, genome, chromosomes);
  std::vector<Pixel<std::int32_t>> buffer(1'000'000);
  buffer.clear();

  auto sel = hf.fetch();

  auto t1 = std::chrono::steady_clock::now();
  std::for_each(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), [&](auto pixel) {
    buffer.emplace_back(std::move(pixel));
    if (buffer.size() == buffer.capacity()) {
      cf.append_pixels(buffer.begin(), buffer.end());
      const auto t2 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()) /
          1000.0;
      spdlog::info(FMT_STRING("Processing {:ucsc} at {:.0f} pixels/s (cache hit rate {:.2f}%)..."),
                   buffer.back().coords.bin1, double(buffer.size()) / delta,
                   hf.block_cache_hit_rate() * 100);
      hf.reset_cache_stats();
      t1 = t2;
      buffer.clear();
    }
  });

  if (!buffer.empty()) {
    cf.append_pixels(buffer.begin(), buffer.end());
  }

  for (const auto norm : normalization_methods) {
    copy_weights(hf, cf, norm, fail_if_norm_not_found);
  }

  t1 = std::chrono::steady_clock::now();
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

  if (c.resolutions.size() == 1) {
    hic::HiCFile hf(c.input_hic.string(), c.resolutions.front(), hic::MatrixType::observed,
                    hic::MatrixUnit::BP);
    convert_resolution(hf, c.output_cooler.string(), chroms, c.resolutions.front(), c.genome,
                       c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable);
    return;
  }

  std::for_each(c.resolutions.rbegin(), c.resolutions.rend(), [&](const auto res) {
    hic::HiCFile hf(c.input_hic.string(), res, hic::MatrixType::observed, hic::MatrixUnit::BP);
    convert_resolution(
        hf, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_cooler.string(), res), chroms,
        res, c.genome, c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable);
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
