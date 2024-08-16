// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./dump.hpp"

#include <fmt/compile.h>
#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/hic.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"

namespace hictk::tools {

template <typename PixelIt>
static void dump_pixels(PixelIt first_pixel, PixelIt last_pixel,
                        const std::shared_ptr<const BinTable>& bins, bool join) {
  if (!join) {
    print_pixels(first_pixel, last_pixel);
    return;
  }
  auto jsel = transformers::JoinGenomicCoords(first_pixel, last_pixel, bins);
  print_pixels(jsel.begin(), jsel.end());
}

static void dump_pixels_gw(File& f, std::string_view normalization, bool join, bool sorted) {
  if (f.is_hic()) {
    f.get<hic::File>().optimize_cache_size_for_iteration();
  }

  if (f.is_hic()) {
    const auto& ff = f.get<hic::File>();
    const auto sel = ff.fetch(balancing::Method{normalization});
    dump_pixels(sel.template begin<double>(sorted), sel.template end<double>(), ff.bins_ptr(),
                join);
    return;
  }

  const auto& ff = f.get<cooler::File>();
  const auto sel = ff.fetch(balancing::Method{normalization});
  dump_pixels(sel.template begin<double>(), sel.template end<double>(), ff.bins_ptr(), join);
}

static void dump_pixels_chrom_chrom(File& f, std::string_view range1, std::string_view range2,
                                    std::string_view normalization, bool join, bool sorted) {
  if (f.is_hic()) {
    const auto& ff = f.get<hic::File>();
    const auto sel = ff.fetch(range1, range2, balancing::Method{normalization});
    dump_pixels(sel.template begin<double>(sorted), sel.template end<double>(), ff.bins_ptr(),
                join);
    return;
  }

  const auto& ff = f.get<cooler::File>();
  const auto sel = ff.fetch(range1, range2, balancing::Method{normalization});
  dump_pixels(sel.template begin<double>(), sel.template end<double>(), ff.bins_ptr(), join);
}

static void dump_pixels(File& f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, bool join, bool sorted) {
  auto norm = balancing::Method{std::string{normalization}};
  if (range1 == "all") {
    assert(range2 == "all");
    dump_pixels_gw(f, normalization, join, sorted);
    return;
  }

  dump_pixels_chrom_chrom(f, range1, range2, normalization, join, sorted);
}

static void process_query_cis_only(File& f, std::string_view normalization, bool join,
                                   bool sorted) {
  for (const auto& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    dump_pixels(f, chrom.name(), chrom.name(), normalization, join, sorted);
  }
}

template <typename File, typename Selector = decltype(std::declval<File>().fetch("chr1", "chr2")),
          typename PixelIt = decltype(std::declval<Selector>().template begin<double>())>
static transformers::PixelMerger<PixelIt> init_pixel_merger(const File& f,
                                                            std::string_view normalization) {
  std::vector<PixelIt> heads{};
  std::vector<PixelIt> tails{};

  for (std::uint32_t chrom1_id = 0; chrom1_id < f.chromosomes().size(); ++chrom1_id) {
    const auto& chrom1 = f.chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < f.chromosomes().size(); ++chrom2_id) {
      const auto& chrom2 = f.chromosomes().at(chrom2_id);
      if (chrom2.is_all()) {
        continue;
      }
      try {
        const auto sel = f.fetch(chrom1.name(), chrom2.name(), balancing::Method{normalization});
        heads.emplace_back(sel.template begin<double>());
        tails.emplace_back(sel.template end<double>());
      } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        const auto missing_norm = msg.find("unable to find") != std::string_view::npos &&
                                  msg.find("normalization vector") != std::string_view::npos;
        if (!missing_norm) {
          throw;
        }
      }
    }
  }

  if (heads.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find {} normalization vectors at {} ({})"), normalization,
                    f.resolution(), hic::MatrixUnit::BP));
  }

  return {std::move(heads), std::move(tails)};
}

static void dump_pixels_trans_only_sorted(File& f, std::string_view normalization, bool join) {
  auto norm = balancing::Method{std::string{normalization}};

  std::visit(
      [&](const auto& ff) {
        const auto merger = init_pixel_merger(ff, normalization);

        if (!join) {
          print_pixels(merger.begin(), merger.end());
        }

        auto jsel = transformers::JoinGenomicCoords(merger.begin(), merger.end(), f.bins_ptr());
        print_pixels(jsel.begin(), jsel.end());
      },
      f.get());
}

static void dump_pixels_trans_only_unsorted(File& f, std::string_view normalization, bool join) {
  const auto& chromosomes = f.chromosomes();
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes.size(); ++chrom1_id) {
    const auto& chrom1 = chromosomes[chrom1_id];
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < chromosomes.size(); ++chrom2_id) {
      const auto& chrom2 = chromosomes[chrom2_id];
      if (chrom2.is_all()) {
        continue;
      }

      dump_pixels(f, chrom1.name(), chrom2.name(), normalization, join, false);
    }
  }
}

static void process_query(File& f, std::string_view table, std::string_view range1,
                          std::string_view range2, std::string_view normalization, bool join,
                          bool sorted) {
  if (table == "bins") {
    dump_bins(f, range1);
    return;
  }

  if (table == "weights") {
    dump_weights(f, range1);
    return;
  }

  assert(table == "pixels");
  dump_pixels(f, range1, range2, normalization, join, sorted);
}

static void process_query_trans_only(File& f, std::string_view normalization, bool join,
                                     bool sorted) {
  if (sorted) {
    dump_pixels_trans_only_sorted(f, normalization, join);
    return;
  }
  dump_pixels_trans_only_unsorted(f, normalization, join);
}

static void dump_tables(const DumpConfig& c) {
  hictk::File f{c.uri, c.resolution, c.matrix_type, c.matrix_unit};

  if (c.query_file.empty() && !c.cis_only && !c.trans_only) {
    process_query(f, c.table, c.range1, c.range2, c.normalization, c.join, c.sorted);
    return;
  }

  if (c.cis_only) {
    assert(c.table == "pixels");
    process_query_cis_only(f, c.normalization, c.join, c.sorted);
    return;
  }

  if (c.trans_only) {
    assert(c.table == "pixels");
    process_query_trans_only(f, c.normalization, c.join, c.sorted);
    return;
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
}

int dump_subcmd(const DumpConfig& c) {
  if (c.table == "bins" || c.table == "pixels" || c.table == "weights") {
    dump_tables(c);
    return 0;
  }

  if (c.table == "chroms") {
    dump_chroms(c.uri, c.format, c.resolution);
    return 0;
  }

  if (c.table == "resolutions") {
    dump_resolutions(c.uri, c.format, c.resolution);
    return 0;
  }

  if (c.table == "normalizations") {
    dump_normalizations(c.uri, c.format, c.resolution);
    return 0;
  }

  assert(c.table == "cells");

  dump_cells(c.uri, c.format);
  return 0;
}
}  // namespace hictk::tools
