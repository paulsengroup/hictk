// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <variant>

#include "./common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers.hpp"

namespace hictk::tools {

template <typename PixelIt>
static void print_pixels(PixelIt first, PixelIt last) {
  std::for_each(first, last, [&](const auto& pixel) { print(pixel); });
}

void dump_pixels(const cooler::File& f, std::string_view range1, std::string_view range2,
                 std::string_view normalization, bool join) {
  auto weight = f.read_weights(normalization);
  if (range1 == "all") {
    assert(range2 == "all");
    auto sel = f.fetch(weight);
    if (!join) {
      return print_pixels(sel.template begin<double>(), sel.template end<double>());
    }

    auto jsel =
        transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
    return print_pixels(jsel.begin(), jsel.end());
  }

  auto sel = f.fetch(range1, range2, weight);
  if (!join) {
    return print_pixels(sel.template begin<double>(), sel.template end<double>());
  }

  auto jsel = transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
  print_pixels(jsel.begin(), jsel.end());
}

void dump_pixels(const hic::HiCFile& f, std::string_view range1, std::string_view range2,
                 std::string_view normalization, bool join) {
  auto norm = hic::ParseNormStr(std::string{normalization});
  if (range1 == "all") {
    assert(range2 == "all");
    auto sel = f.fetch(norm);
    if (!join) {
      return print_pixels(sel.template begin<double>(), sel.template end<double>());
    }

    auto jsel =
        transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
    return print_pixels(jsel.begin(), jsel.end());
  }

  auto sel = f.fetch(range1, range2, norm);
  if (!join) {
    return print_pixels(sel.template begin<double>(), sel.template end<double>());
  }

  auto jsel = transformers::JoinGenomicCoords(sel.begin<double>(), sel.end<double>(), f.bins_ptr());
  print_pixels(jsel.begin(), jsel.end());
}

template <typename File>
static void process_query(const File& f, std::string_view table, std::string_view range1,
                          std::string_view range2, std::string_view normalization, bool join) {
  if (table == "chroms") {
    return dump_chroms(f, range1);
  }
  if (table == "bins") {
    return dump_bins(f, range1);
  }

  assert(table == "pixels");
  return dump_pixels(f, range1, range2, normalization, join);
}

using FileVar = std::variant<cooler::File, hic::HiCFile>;

[[nodiscard]] static FileVar open_hic_file(std::string_view path, std::uint32_t resolution,
                                           hic::MatrixType matrix_type,
                                           hic::MatrixUnit matrix_unit) {
  return {hic::HiCFile(std::string{path}, resolution, matrix_type, matrix_unit)};
}

[[nodiscard]] static FileVar open_cooler_file(std::string_view uri) {
  return {cooler::File::open_read_only_read_once(uri)};
}

void dump_subcmd(const DumpConfig& c) {
  const auto is_cooler = cooler::utils::is_cooler(c.uri) || cooler::utils::is_multires_file(c.uri);
  auto file{is_cooler ? open_cooler_file(c.uri)
                      : open_hic_file(c.uri, c.resolution, c.matrix_type, c.matrix_unit)};

  std::visit(
      [&](const auto& f) {
        if (c.query_file.empty()) {
          process_query(f, c.table, c.range1, c.range2, c.normalization, c.join);
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
          process_query(f, c.table, range1, range2, c.normalization, c.join);
        }
      },
      file);
}
}  // namespace hictk::tools
