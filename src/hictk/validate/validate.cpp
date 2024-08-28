// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

#include "hictk/chromosome.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

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
        fmt::format(FMT_STRING("### FAILURE: \"{}\" is not a valid .hic file:\n"
                               "Validation failed for {}:{} map at {} resolution:\n"
                               "{}"),
                    hf.path(), chrom1.name(), chrom2.name(), hf.resolution(), e.what()));
  }
}

static int validate_hic(const std::string& path, bool quiet) {
  try {
    for (const auto& res : hic::utils::list_resolutions(path)) {
      const hic::File hf(path, res);
      const auto& chroms = hf.chromosomes();
      for (std::uint32_t i = 0; i < chroms.size(); ++i) {
        for (std::uint32_t j = i; j < chroms.size(); ++j) {
          const auto& chrom1 = chroms.at(i);
          const auto& chrom2 = chroms.at(j);

          validate_hic(hf, chrom1, chrom2);
        }
      }
    }
  } catch (const std::exception& e) {
    if (!quiet) {
      fmt::print(FMT_STRING("{}\n"), e.what());
    }
    return 1;
  }

  if (!quiet) {
    fmt::print(FMT_STRING("### SUCCESS: \"{}\" is a valid .hic file.\n"), path);
  }
  return 0;
}

[[nodiscard]] static bool validate_bin_table_size(const cooler::File& clr, bool quiet) {
  const auto& chroms = clr.dataset("bins/chrom");
  const auto& starts = clr.dataset("bins/start");
  const auto& ends = clr.dataset("bins/end");

  const auto expected_num_bins = clr.bins().size();
  if (chroms.size() != expected_num_bins || starts.size() != expected_num_bins ||
      ends.size() != expected_num_bins) {
    if (!quiet) {
      fmt::print(FMT_STRING("bin_table_has_correct_length=false\n"));
      return false;
    }
  }
  if (!quiet) {
    fmt::print(FMT_STRING("bin_table_has_correct_length=true\n"));
  }

  return true;
}

[[nodiscard]] static bool validate_bins_dtypes(const cooler::File& clr, bool quiet) {
  try {
    std::ignore = *clr.dataset("bins/chrom").begin<std::string>();
    std::ignore = *clr.dataset("bins/start").begin<std::int32_t>();
    std::ignore = *clr.dataset("bins/end").begin<std::int32_t>();
    if (!quiet) {
      fmt::print(FMT_STRING("bin_table_has_correct_dtypes=true\n"));
    }
    return true;
  } catch (...) {
    if (!quiet) {
      fmt::print(FMT_STRING("bin_table_has_correct_dtypes=false\n"));
    }
    return false;
  }
}

[[nodiscard]] static bool validate_bins(const cooler::File& clr, bool quiet) {
  const auto& chroms = clr.dataset("bins/chrom");
  const auto& starts = clr.dataset("bins/start");
  const auto& ends = clr.dataset("bins/end");

  auto first_chrom_id = chroms.begin<std::int32_t>();
  const auto last_chrom = chroms.end<std::int32_t>();
  auto first_start = starts.begin<std::int32_t>();
  const auto last_start = starts.end<std::int32_t>();
  auto first_end = ends.begin<std::int32_t>();
  const auto last_end = ends.end<std::int32_t>();

  auto first_bin = clr.bins().begin();
  const auto last_bin = clr.bins().end();

  std::size_t num_invalid_bins{};
  const auto num_bins = clr.bins().size();
  for (std::size_t i = 0; i < num_bins; ++i) {
    const auto bin = *first_bin++;

    const auto chrom_it = clr.chromosomes().find(static_cast<std::uint32_t>(*first_chrom_id++));
    if (chrom_it == clr.chromosomes().end()) {
      ++first_start;
      ++first_end;
      ++num_invalid_bins;
      continue;
    }

    const auto& chrom = *chrom_it;
    const auto start = static_cast<std::uint32_t>(*first_start++);
    const auto end = static_cast<std::uint32_t>(*first_end++);

    if (bin.chrom() != chrom || bin.start() != start || bin.end() != end) {
      ++num_invalid_bins;
    }
  }

  if (!quiet) {
    fmt::print(FMT_STRING("num_invalid_bins={}"), num_invalid_bins);
  }

  return num_invalid_bins == 0;
}

[[nodiscard]] static bool check_bin_table(const cooler::File& clr, bool quiet) {
  if (!validate_bin_table_size(clr, quiet)) {
    return false;
  }
  if (!validate_bins_dtypes(clr, quiet)) {
    return false;
  }

  return validate_bins(clr, quiet);
}

static int validate_cooler(std::string_view path, bool validate_index, bool quiet) {
  auto status = cooler::utils::is_cooler(path);
  std::optional<cooler::File> clr{};
  if (status.is_cooler) {
    try {
      clr = cooler::File(path);
    } catch (...) {
      status.is_cooler = false;
    }
  }

  std::optional<bool> bins_ok{};
  if (status.is_cooler) {
    assert(clr.has_value());
    bins_ok = check_bin_table(*clr, quiet);
  }

  auto index_ok = true;
  if (status.is_cooler && validate_index) {
    index_ok = cooler::utils::index_is_valid(path, quiet);
  }

  const auto cooler_is_valid = !!status && (bins_ok.has_value() && *bins_ok) && index_ok;

  if (!quiet) {
    fmt::print(FMT_STRING("{}\n"), status);
    if (validate_index) {
      fmt::print(FMT_STRING("index_is_valid={}\n"), index_ok);
    } else {
      fmt::print(FMT_STRING("index_is_valid=not_checked\n"));
    }
    if (bins_ok.has_value()) {
      fmt::print(FMT_STRING("bin_table_is_valid={}\n"), *bins_ok);
    } else {
      fmt::print(FMT_STRING("bin_table_is_valid=not_checked\n"));
    }

    fmt::print(FMT_STRING("### {}: \"{}\" {} a valid Cooler.\n"),
               cooler_is_valid ? "SUCCESS" : "FAILURE", path, cooler_is_valid ? "is" : "is not");
  }
  return static_cast<int>(!cooler_is_valid);
}

static int validate_mcool(std::string_view path, bool validate_index, bool quiet) {
  const cooler::MultiResFile mclr{std::filesystem::path(path)};
  auto resolutions = mclr.resolutions();
  std::sort(resolutions.begin(), resolutions.end(), std::greater{});
  for (const auto& res : resolutions) {
    const auto status = validate_cooler(mclr.open(res).uri(), validate_index, quiet);
    if (status != 0) {
      return 1;
    }
  }
  return 0;
}

static int validate_scool(std::string_view path, bool validate_index, bool quiet) {
  const cooler::SingleCellFile sclr{std::filesystem::path(path)};
  for (const auto& cell : sclr.cells()) {
    const auto status = validate_cooler(sclr.open(cell).uri(), validate_index, quiet);
    if (status != 0) {
      return 1;
    }
  }
  return 0;
}

int validate_subcmd(const ValidateConfig& c) {
  try {
    const auto is_cooler = cooler::utils::is_cooler(c.uri);
    const auto is_hic = hic::utils::is_hic_file(c.uri);
    const auto is_mcool = cooler::utils::is_multires_file(c.uri);
    const auto is_scool = cooler::utils::is_scool_file(c.uri);

    if (!is_hic && !is_cooler && !is_mcool && !is_scool) {
      if (!c.quiet) {
        fmt::print(FMT_STRING("### FAILURE: \"{}\" is not in .hic or .[ms]cool format!\n"), c.uri);
      }
      return 1;
    }

    if (is_hic) {
      return validate_hic(c.uri, c.quiet);
    }

    if (is_mcool) {
      return validate_mcool(c.uri, c.validate_index, c.quiet);
    }

    if (is_scool) {
      return validate_scool(c.uri, c.validate_index, c.quiet);
    }

    return validate_cooler(c.uri, c.validate_index, c.quiet);
  } catch (...) {
    if (c.quiet) {
      return 1;
    }
    throw;
  }
}

}  // namespace hictk::tools
