// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <variant>

#include "./common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

template <typename N, bool join>
static void print_pixels(typename cooler::PixelSelector<N>::iterator first_pixel,
                         typename cooler::PixelSelector<N>::iterator last_pixel,
                         const std::shared_ptr<const cooler::Weights>& weights) {
  if (weights != nullptr) {
    const auto sel = cooler::Balancer<N>(first_pixel, last_pixel, weights);
    std::for_each(sel.begin(), sel.end(), [&](const auto& pixel) { print<join>(pixel); });
    return;
  }
  std::for_each(first_pixel, last_pixel, [&](const auto& pixel) { print<join>(pixel); });
}

template <bool join, typename HiCPixelIt>
static void print_pixels(HiCPixelIt first_pixel, HiCPixelIt last_pixel) {
  std::for_each(first_pixel, last_pixel, [&](const auto& pixel) { print<join>(pixel); });
}

template <bool join>
static void dump_pixels(const cooler::File& clr, std::string_view range1, std::string_view range2,
                        std::string_view normalization) {
  const auto has_int_pixels = clr.has_integral_pixels();
  const auto weights =  // TODO: pass ptr to dump_pixels
      normalization == "NONE" ? std::shared_ptr<const cooler::Weights>(nullptr)
                              : clr.read_weights(normalization);

  if (range1 == "all") {
    assert(range2 == "all");
    return has_int_pixels
               ? print_pixels<std::int64_t, join>(clr.begin<std::int64_t>(),
                                                  clr.end<std::int64_t>(), weights)
               : print_pixels<double, join>(clr.begin<double>(), clr.end<double>(), weights);
  }

  if (has_int_pixels) {
    auto sel = clr.fetch<std::int64_t>(range1, range2);
    return print_pixels<std::int64_t, join>(sel.begin(), sel.end(), weights);
  }

  auto sel = clr.fetch<double>(range1, range2);
  return print_pixels<double, join>(sel.begin(), sel.end(), weights);
}

template <bool join>
static void dump_pixels(const hic::HiCFile& f, std::string_view range1,
                        [[maybe_unused]] std::string_view range2, std::string_view normalization) {
  if (range1 == "all") {
    assert(range2 == "all");
    auto sel = f.fetch(hic::ParseNormStr(std::string{normalization}));
    return print_pixels<join>(sel.begin<double>(), sel.end<double>());
  }
  auto sel = f.fetch(range1, range2);
  return print_pixels<join>(sel.begin<double>(), sel.end<double>());
}

template <bool join, typename File>
static void process_query(const File& f, std::string_view table, std::string_view range1,
                          std::string_view range2, std::string_view normalization) {
  if (table == "chroms") {
    return dump_chroms(f, range1);
  }
  if (table == "bins") {
    return dump_bins(f, range1);
  }

  assert(table == "pixels");
  return dump_pixels<join>(f, range1, range2, normalization);
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
          c.join ? process_query<true>(f, c.table, c.range1, c.range2, c.normalization)
                 : process_query<false>(f, c.table, c.range1, c.range2, c.normalization);
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
          c.join ? process_query<true>(f, c.table, range1, range2, c.normalization)
                 : process_query<false>(f, c.table, range1, range2, c.normalization);
        }
      },
      file);
}
}  // namespace hictk::tools
