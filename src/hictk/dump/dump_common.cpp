// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "./dump.hpp"
#include "hictk/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {
void print(const Pixel<double>& pixel) {
  fmt::print(FMT_COMPILE("{:bg2}\t{:.16g}\n"), pixel.coords, pixel.count);
}
void print(const ThinPixel<double>& pixel) {
  fmt::print(FMT_COMPILE("{:d}\t{:d}\t{:.16g}\n"), pixel.bin1_id, pixel.bin2_id, pixel.count);
}

void dump_bins(const File& f, std::string_view range1, std::string_view range2) {
  if (range1 == "all") {
    assert(range2 == "all");
    for (const auto& bin : f.bins()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
    }
    return;
  }

  auto coords1 = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range1});
  auto coords2 = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range2});
  if (coords1 > coords2) {
    std::swap(coords1, coords2);
  }
  auto [first_bin1, last_bin1] = f.bins().find_overlap(coords1);
  std::for_each(first_bin1, last_bin1, [](const Bin& bin) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
  });

  if (coords1 == coords2) {
    return;
  }

  auto [first_bin2, last_bin2] = f.bins().find_overlap(coords2);
  std::for_each(first_bin2, last_bin2, [](const Bin& bin) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
  });
}

[[nodiscard]] static std::pair<std::size_t, std::size_t> compute_bin_ids(const BinTable& bins,
                                                                         std::string_view range) {
  if (range == "all") {
    return {0, bins.size()};
  }

  const auto coords = GenomicInterval::parse_ucsc(bins.chromosomes(), std::string{range});
  const auto [first_bin, last_bin] = bins.find_overlap(coords);
  const auto i0 = (*first_bin).id();
  const auto i1 = (*last_bin).id();

  return {i0, i1};
}

static void dump_weights(const BinTable& bins, std::string_view range,
                         const std::vector<balancing::Method>& norms,
                         const std::vector<balancing::Weights>& weights, bool print_header) {
  const auto [i0, i1] = compute_bin_ids(bins, range);

  if (print_header) {
    fmt::print(FMT_STRING("{}\n"), fmt::join(norms, "\t"));
  }
  std::vector<double> record(norms.size());
  for (std::size_t i = i0; i < i1; ++i) {
    for (std::size_t j = 0; j < norms.size(); ++j) {
      const auto& w = weights[j];
      if (w.type() == balancing::Weights::Type::DIVISIVE) {
        record[j] = w[i];
      } else {
        record[j] = 1.0 / w[i];
      }
    }
    fmt::print(FMT_COMPILE("{}\n"), fmt::join(record, "\t"));
  }
}

void dump_weights(const File& f, std::string_view range1, std::string_view range2) {
  const auto norms = f.avail_normalizations();
  if (norms.empty()) {
    return;
  }

  std::vector<balancing::Weights> weights{};
  for (const auto& norm : norms) {
    weights.emplace_back(f.normalization(norm.to_string()));  // NOLINT
  }

  if (range1 == "all") {
    assert(range2 == "all");
    dump_weights(f.bins(), range1, norms, weights, true);
    return;
  }

  auto coords1 = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range1});
  auto coords2 = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range2});
  if (coords1 > coords2) {
    std::swap(range1, range2);
    std::swap(coords1, coords2);
  }

  dump_weights(f.bins(), range1, norms, weights, true);

  if (range1 != range2) {
    dump_weights(f.bins(), range2, norms, weights, false);
  }
}

void dump_cells(std::string_view uri, std::string_view format) {
  if (format != "scool") {
    throw std::runtime_error(fmt::format(FMT_STRING("\"{}\" is not a .scool file"), uri));
  }
  const auto cells = cooler::SingleCellFile{uri}.cells();
  std::for_each(cells.begin(), cells.end(),
                [&](const auto& cell) { fmt::print(FMT_COMPILE("{}\n"), cell); });
}

void dump_chroms(std::string_view uri, std::string_view range1, std::string_view range2,
                 std::string_view format, std::optional<std::uint32_t> resolution) {
  Reference ref{};

  if (format == "mcool") {
    ref = cooler::MultiResFile{std::string{uri}}.chromosomes();
  } else if (format == "scool") {
    ref = cooler::SingleCellFile{std::string{uri}}.chromosomes();
  } else {
    ref = File{std::string{uri}, resolution}.chromosomes();
  }

  if (range1 == "all") {
    assert(range2 == "all");

    for (const Chromosome& chrom : ref) {
      if (!chrom.is_all()) {
        fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), chrom.name(), chrom.size());
      }
    }
    return;
  }

  auto coords1 = GenomicInterval::parse_ucsc(ref, std::string{range1});
  auto coords2 = GenomicInterval::parse_ucsc(ref, std::string{range2});

  if (coords1 > coords2) {
    std::swap(coords1, coords2);
  }

  fmt::print(FMT_STRING("{:s}\t{:d}\n"), coords1.chrom().name(), coords1.chrom().size());
  if (coords1 == coords2) {
    return;
  }
  fmt::print(FMT_STRING("{:s}\t{:d}\n"), coords2.chrom().name(), coords2.chrom().size());
}

static phmap::btree_set<std::string> get_normalizations(std::string_view uri,
                                                        std::string_view format,
                                                        std::optional<std::uint32_t> resolution) {
  assert(format != "mcool");
  if (format == "scool") {
    const auto cell_ids = cooler::SingleCellFile{uri}.cells();
    if (cell_ids.empty()) {
      return {};
    }

    const auto scool_uri = fmt::format(FMT_STRING("{}::/cells/{}"), uri, *cell_ids.begin());
    return get_normalizations(scool_uri, "cool", std::nullopt);
  }

  phmap::btree_set<std::string> norms{};
  const auto norms_ = File{std::string{uri}, resolution}.avail_normalizations();
  std::transform(norms_.begin(), norms_.end(), std::inserter(norms, norms.begin()),
                 [](const auto& n) { return std::string{n.to_string()}; });

  return norms;
}

void dump_normalizations(std::string_view uri, std::string_view format,
                         std::optional<std::uint32_t> resolution) {
  phmap::btree_set<std::string> norms{};
  std::vector<std::uint32_t> resolutions{};
  if (format == "mcool") {
    resolutions = cooler::MultiResFile{uri}.resolutions();
    if (resolutions.empty()) {
      return;
    }
  } else if (format == "hic" && !resolution.has_value()) {
    resolutions = hic::utils::list_resolutions(std::string{uri});
    if (resolutions.empty()) {
      return;
    }
  }

  if (resolutions.empty()) {
    norms = get_normalizations(uri, format, resolution);
  } else {
    format = format == "hic" ? "hic" : "cool";
    std::for_each(resolutions.begin(), resolutions.end(),
                  [&](const auto res) { norms.merge(get_normalizations(uri, format, res)); });
  }

  if (!norms.empty()) {
    fmt::print(FMT_STRING("{}\n"), fmt::join(norms, "\n"));
  }
}

void dump_resolutions(std::string_view uri, std::string_view format,
                      std::optional<std::uint32_t> resolution) {
  std::vector<std::uint32_t> resolutions{};

  if (format == "hic") {
    resolutions = hic::utils::list_resolutions(uri);
    if (resolution.has_value()) {
      const auto res_found =
          std::find(resolutions.begin(), resolutions.end(), *resolution) != resolutions.end();
      resolutions.clear();
      if (res_found) {
        resolutions.push_back(*resolution);
      }
    }
  } else if (format == "mcool") {
    resolutions = cooler::MultiResFile{uri}.resolutions();
  } else if (format == "scool") {
    resolutions.push_back(cooler::SingleCellFile{uri}.resolution());
  } else {
    assert(format == "cool");
    resolutions.push_back(cooler::File{uri}.resolution());
  }

  if (!resolutions.empty()) {
    fmt::print(FMT_STRING("{}\n"), fmt::join(resolutions, "\n"));
  }
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
std::pair<std::string, std::string> parse_bedpe(std::string_view line) {
  if (line.empty()) {
    throw std::runtime_error("found an empty line");
  }

  if (line.back() == '\r') {
    line = line.substr(0, line.size() - 1);
  }

  auto next_token = [&]() {
    assert(!line.empty());
    const auto pos1 = line.find('\t');
    const auto pos2 = line.find('\t', pos1 + 1);
    const auto pos3 = line.find('\t', pos2 + 1);

    auto tok = std::string{line.substr(0, pos3)};
    tok[pos1] = ':';
    tok[pos2] = '-';
    line.remove_prefix(pos3 + 1);
    return tok;
  };
  const auto range1 = next_token();
  const auto range2 = next_token();

  return std::make_pair(range1, range2);
}
}  // namespace hictk::tools
