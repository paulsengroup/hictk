// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <variant>

#include "hictk/balancing/methods.hpp"
#include "hictk/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers.hpp"

namespace hictk::tools {

static void print(const Pixel<double>& pixel) {
  fmt::print(FMT_COMPILE("{:bg2}\t{:.16g}\n"), pixel.coords, pixel.count);
}
static void print(const ThinPixel<double>& pixel) {
  fmt::print(FMT_COMPILE("{:d}\t{:d}\t{:.16g}\n"), pixel.bin1_id, pixel.bin2_id, pixel.count);
}

template <typename File>
static void dump_bins(const File& f, std::string_view range) {
  if (range == "all") {
    for (const auto& bin : f.bins()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
    }
    return;
  }

  const auto coords = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range});
  auto [first_bin, last_bin] = f.bins().find_overlap(coords);
  std::for_each(first_bin, last_bin, [](const Bin& bin) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
  });
}

[[nodiscard]] static std::pair<std::string, std::string> parse_bedpe(std::string_view line) {
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

template <typename PixelIt>
static void print_pixels(PixelIt first, PixelIt last) {
  std::for_each(first, last, [&](const auto& pixel) { print(pixel); });
}

static void dump_pixels(const cooler::File& f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, bool join, [[maybe_unused]] bool sorted) {
  auto weights = f.read_weights(balancing::Method{normalization});
  if (range1 == "all") {
    assert(range2 == "all");
    auto sel = f.fetch(weights);
    if (!join) {
      return print_pixels(sel.template begin<double>(), sel.template end<double>());
    }

    auto jsel =
        transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
    return print_pixels(jsel.begin(), jsel.end());
  }

  auto sel = f.fetch(range1, range2, weights);
  if (!join) {
    return print_pixels(sel.template begin<double>(), sel.template end<double>());
  }

  auto jsel = transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
  print_pixels(jsel.begin(), jsel.end());
}

static void dump_pixels(hic::File& f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, bool join, bool sorted) {
  auto norm = balancing::Method{std::string{normalization}};
  if (range1 == "all") {
    assert(range2 == "all");
    f.optimize_cache_size_for_iteration();
    auto sel = f.fetch(norm);
    if (!join) {
      return print_pixels(sel.template begin<double>(sorted), sel.template end<double>());
    }
    auto jsel =
        transformers::JoinGenomicCoords(sel.begin<double>(sorted), sel.end<double>(), f.bins_ptr());
    return print_pixels(jsel.begin(), jsel.end());
  }

  auto sel = f.fetch(range1, range2, norm);
  if (!join) {
    return print_pixels(sel.template begin<double>(sorted), sel.template end<double>());
  }

  auto jsel =
      transformers::JoinGenomicCoords(sel.begin<double>(sorted), sel.end<double>(), f.bins_ptr());
  print_pixels(jsel.begin(), jsel.end());
}

static void process_query(File& f, std::string_view table, std::string_view range1,
                          std::string_view range2, std::string_view normalization, bool join,
                          bool sorted) {
  if (table == "bins") {
    return dump_bins(f, range1);
  }

  assert(table == "pixels");
  std::visit([&](auto& ff) { return dump_pixels(ff, range1, range2, normalization, join, sorted); },
             f.get());
}

static int dump_chroms(std::string_view uri, std::string_view format, std::uint32_t resolution) {
  Reference ref{};

  if (format == "mcool") {
    ref = cooler::MultiResFile{std::string{uri}}.chromosomes();
  } else if (format == "scool") {
    ref = cooler::SingleCellFile{std::string{uri}}.chromosomes();
  } else {
    ref = File{std::string{uri}, resolution}.chromosomes();
  }

  for (const Chromosome& chrom : ref) {
    if (!chrom.is_all()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), chrom.name(), chrom.size());
    }
  }
  return 0;
}

static phmap::btree_set<std::string> get_normalizations(std::string_view uri,
                                                        std::string_view format,
                                                        std::uint32_t resolution) {
  assert(format != "mcool");
  assert(format != "hic" || resolution != 0);
  if (format == "scool") {
    const auto cell_ids = cooler::SingleCellFile{uri}.cells();
    if (cell_ids.empty()) {
      return {};
    }

    const auto scool_uri = fmt::format(FMT_STRING("{}::/cells/{}"), uri, *cell_ids.begin());
    return get_normalizations(scool_uri, "cool", 0);
  }

  phmap::btree_set<std::string> norms{};
  if (uri == "hic" && resolution == 0) {
    const hic::File hf{std::string{uri}, resolution};

    for (const auto& norm : hf.avail_normalizations()) {
      norms.emplace(std::string{norm.to_string()});
    }
    return norms;
  }

  const auto norms_ = File{std::string{uri}, resolution}.avail_normalizations();
  std::transform(norms_.begin(), norms_.end(), std::inserter(norms, norms.begin()),
                 [](const auto& n) { return std::string{n.to_string()}; });

  return norms;
}

static int dump_normalizations(std::string_view uri, std::string_view format,
                               std::uint32_t resolution) {
  phmap::btree_set<std::string> norms{};
  std::vector<std::uint32_t> resolutions{};
  if (format == "mcool") {
    resolutions = cooler::MultiResFile{uri}.resolutions();
    if (resolutions.empty()) {
      return 0;
    }
  } else if (format == "hic" && resolution == 0) {
    resolutions = hic::utils::list_resolutions(std::string{uri});
    if (resolutions.empty()) {
      return 0;
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
  return 0;
}

static int dump_resolutions(std::string_view uri, std::string_view format,
                            std::uint32_t resolution) {
  std::vector<std::uint32_t> resolutions{};

  if (format == "hic") {
    resolutions = hic::utils::list_resolutions(uri);
    if (resolution != 0) {
      const auto res_found =
          std::find(resolutions.begin(), resolutions.end(), resolution) != resolutions.end();
      resolutions.clear();
      if (res_found) {
        resolutions.push_back(resolution);
      }
    }
  } else if (format == "mcool") {
    resolutions = cooler::MultiResFile{uri}.resolutions();
  } else if (format == "scool") {
    resolutions.push_back(cooler::SingleCellFile{uri}.bin_size());
  } else {
    assert(format == "cool");
    resolutions.push_back(cooler::File{uri}.bin_size());
  }

  if (!resolutions.empty()) {
    fmt::print(FMT_STRING("{}\n"), fmt::join(resolutions, "\n"));
  }
  return 0;
}

static int dump_cells(std::string_view uri, std::string_view format) {
  if (format != "scool") {
    throw std::runtime_error(fmt::format(FMT_STRING("\"{}\" is not a .scool file"), uri));
  }
  const auto cells = cooler::SingleCellFile{uri}.cells();
  if (!cells.empty()) {
    fmt::print(FMT_STRING("{}\n"), fmt::join(cells, "\n"));
  }
  return 0;
}

static int dump_tables(const DumpConfig& c) {
  hictk::File f{c.uri, c.resolution, c.matrix_type, c.matrix_unit};

  if (c.query_file.empty()) {
    process_query(f, c.table, c.range1, c.range2, c.normalization, c.join, c.sorted);
    return 0;
  }

  const auto read_from_stdin = c.query_file == "-";
  std::ifstream ifs{};
  ifs.exceptions(ifs.exceptions() | std::ios_base::badbit | std::ios_base::failbit);

  if (!read_from_stdin) {
    assert(std::filesystem::exists(c.query_file));
    ifs.open(c.query_file);
  }

  std::string line;
  while (std::getline(read_from_stdin ? std::cin : ifs, line)) {
    const auto [range1, range2] = parse_bedpe(line);
    process_query(f, c.table, range1, range2, c.normalization, c.join, c.sorted);
  }

  return 0;
}

int dump_subcmd(const DumpConfig& c) {
  if (c.table == "bins" || c.table == "pixels") {
    return dump_tables(c);
  }

  if (c.table == "chroms") {
    return dump_chroms(c.uri, c.format, c.resolution);
  }

  if (c.table == "resolutions") {
    return dump_resolutions(c.uri, c.format, c.resolution);
  }

  if (c.table == "normalizations") {
    return dump_normalizations(c.uri, c.format, c.resolution);
  }

  assert(c.table == "cells");

  return dump_cells(c.uri, c.format);
}
}  // namespace hictk::tools
