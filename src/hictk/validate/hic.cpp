// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "./validate.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

static std::optional<hic::File> open_hic_noexcept(const std::string& path,
                                                  std::uint32_t resolution) noexcept {
  try {
    return hic::File(path, resolution);
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("[{}] failed to open file: {}"), resolution, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("[{}] failed to open file: unknown error"), resolution);
    return {};
  }
}

static void validate_hic(const hic::File& hf, const Chromosome& chrom1, const Chromosome& chrom2) {
  if (chrom1.is_all() || chrom2.is_all()) {
    return;
  }

  try {
    std::ignore = hf.fetch(chrom1.name(), chrom2.name());
  } catch (const std::exception& e) {
    const std::string_view msg{e.what()};
    if (msg.find("Unable to find block map") != std::string_view::npos) {
      return;
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Validation failed for {}:{} map at {} resolution: {}"),
                    chrom1.name(), chrom2.name(), hf.resolution(), e.what()));
  }
}

std::pair<int, toml::table> validate_hic(const std::string& path, bool exhaustive) {
  toml::table status;
  int return_code = 0;
  for (const auto& res : hic::utils::list_resolutions(path)) {
    const auto hf = open_hic_noexcept(path, res);
    if (!hf) {
      status.insert(fmt::to_string(res), "unable to open resolution");
      status.insert_or_assign("is_valid_hic", false);
      return_code = 1;
      // assert(false);  // This should never happen
      if (!exhaustive) {
        return std::make_pair(return_code, status);
      }
      continue;
    }

    const auto& chroms = hf->chromosomes();
    for (std::uint32_t i = 0; i < chroms.size(); ++i) {
      for (std::uint32_t j = i; j < chroms.size(); ++j) {
        const auto& chrom1 = chroms.at(i);
        const auto& chrom2 = chroms.at(j);

        try {
          validate_hic(*hf, chrom1, chrom2);  // NOLINT(bugprone-unchecked-optional-access)
        } catch ([[maybe_unused]] const std::exception& e) {
          SPDLOG_DEBUG(FMT_STRING("[{}]: validation failed for {}:{} {}"), res, chrom1.name(),
                       chrom2.name(), e.what());
          status.insert(fmt::format(FMT_STRING("{}:{}_{}"), chrom1.name(), chrom2.name(), res),
                        "unable to fetch interactions");
          status.insert_or_assign("is_valid_hic", false);
          return_code = 1;
          if (!exhaustive) {
            return std::make_pair(return_code, status);
          }
        }
      }
    }
  }

  status.insert_or_assign("is_valid_hic", return_code == 0);

  return std::make_pair(return_code, status);
}

}  // namespace hictk::tools
