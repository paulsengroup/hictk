// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cassert>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void copy_pixels(hic::internal::HiCFileWriter& w, const cooler::File& base_clr,
                        const ConvertConfig& c) {
  if (c.input_format == "cool") {
    w.add_pixels(base_clr.bin_size(), base_clr.begin<float>(), base_clr.end<float>());
    return;
  }

  assert(c.input_format == "mcool");
  const cooler::MultiResFile mclr(c.path_to_input.string());

  for (const auto& res : c.resolutions) {
    try {
      const auto clr = mclr.open(res);
      w.add_pixels(res, clr.begin<float>(), clr.end<float>());
    } catch (const std::exception& e) {
      const std::string_view msg{e.what()};
      const auto pos = msg.find("does not have interactions for resolution");
      if (pos == std::string_view::npos) {
        throw;
      }
    }
  }
}

static void copy_normalization_vector(hic::internal::HiCFileWriter& w, const cooler::File& clr,
                                      const balancing::Method& norm, bool throw_if_missing) {
  if (norm == balancing::Method::NONE()) {
    return;
  }

  try {
    const auto& weights = *clr.read_weights(norm);
    std::vector<float> weights_f(weights().size());
    std::transform(weights().begin(), weights().end(), weights_f.begin(), [&](const double n) {
      if (weights.type() == balancing::Weights::Type::MULTIPLICATIVE) {
        return static_cast<float>(1.0 / n);
      }
      return static_cast<float>(n);
    });

    const auto norm_name = norm.to_string() == "weight" ? "ICE" : norm.to_string();
    SPDLOG_INFO(FMT_STRING("[{}] adding {} normalization vector"), clr.bin_size(), norm_name);
    w.add_norm_vector(norm_name, "BP", clr.bin_size(), weights_f, true);

  } catch (const std::exception& e) {
    const std::string_view msg{e.what()};
    const auto match = msg.find(fmt::format(FMT_STRING("unable to read \"{}\" weights"), norm));
    if (match == std::string_view::npos) {
      throw;
    }
    if (throw_if_missing) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find {} normalization vector for resolution {}"), norm,
                      norm, clr.bin_size()));
    }
    SPDLOG_WARN(FMT_STRING("[{}] {} normalization vector is missing. SKIPPING!"), clr.bin_size(),
                norm);
  }
}

static void copy_normalization_vectors(hic::internal::HiCFileWriter& w,
                                       const cooler::File& base_clr, const ConvertConfig& c) {
  const auto avail_normalizations = base_clr.avail_normalizations();

  if (c.input_format == "cool") {
    for (const auto& norm : c.normalization_methods) {
      copy_normalization_vector(w, base_clr, norm, c.fail_if_normalization_method_is_not_avaliable);
    }
    w.write_norm_vectors_and_norm_expected_values();
    return;
  }

  assert(c.input_format == "mcool");
  const cooler::MultiResFile mclr(c.path_to_input.string());

  for (const auto& res : c.resolutions) {
    const auto clr = mclr.open(res);
    for (const auto& norm : c.normalization_methods) {
      copy_normalization_vector(w, clr, norm, c.fail_if_normalization_method_is_not_avaliable);
    }
  }
  w.write_norm_vectors_and_norm_expected_values();
}

void cool_to_hic(const ConvertConfig& c) {
  if (c.force && std::filesystem::exists(c.path_to_output)) {
    [[maybe_unused]] std::error_code ec{};
    std::filesystem::remove(c.path_to_output, ec);
  }

  const auto base_uri = c.input_format == "cool"
                            ? c.path_to_input.string()
                            : fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                          c.path_to_input.string(), c.resolutions.front());

  const cooler::File base_clr{base_uri};
  if (base_clr.bin_size() == 0) {
    throw std::runtime_error("converting cooler files with variable bin size is not supported.");
  }

  const auto chromosomes = cooler::File(base_uri).chromosomes();
  const auto& resolutions = c.resolutions;

  const internal::TmpDir tmpdir{c.tmp_dir};
  hictk::hic::internal::HiCFileWriter w(c.path_to_output.string(), chromosomes, resolutions,
                                        c.genome, c.threads, c.chunk_size, c.tmp_dir,
                                        c.compression_lvl, c.skip_all_vs_all_matrix);
  copy_pixels(w, base_clr, c);
  w.serialize();

  copy_normalization_vectors(w, base_clr, c);
}
}  // namespace hictk::tools
