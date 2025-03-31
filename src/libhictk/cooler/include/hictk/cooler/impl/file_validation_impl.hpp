// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <highfive/H5Exception.hpp>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/fmt/pixel.hpp"
#include "hictk/pixel.hpp"

namespace hictk::cooler {

inline void File::validate_bins(bool full) const {
  try {
    auto nchroms = dataset("bins/chrom").size();
    auto nstarts = dataset("bins/start").size();
    auto nends = dataset("bins/end").size();
    if (nchroms != nstarts || nchroms != nends) {
      throw std::runtime_error(fmt::format(FMT_STRING("Datasets have inconsistent sizes:\n"
                                                      " - \"bins/chrom\": {}\n"
                                                      " - \"bins/start\": {}\n"
                                                      " - \"bins/end\": {}\n"
                                                      "Expected {}"),
                                           nchroms, nstarts, nends, bins().size()));
    }

    const auto &nbins = nchroms;
    if (nbins != bins().size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} bins, found {}"), bins().size(), nchroms));
    }

    // NOLINTBEGIN(*-avoid-magic-numbers)
    auto chrom_it = dataset("bins/chrom").begin<std::uint32_t>(64'000);
    auto start_it = dataset("bins/start").begin<std::uint32_t>(64'000);
    auto end_it = dataset("bins/end").begin<std::uint32_t>(64'000);
    // NOLINTEND(*-avoid-magic-numbers)

    auto last_chrom = dataset("bins/chrom").end<std::uint32_t>(0);
    auto last_start = dataset("bins/start").end<std::uint32_t>(0);
    auto last_end = dataset("bins/end").end<std::uint32_t>(0);

    if (full) {
      std::size_t i = 0;
      for (const Bin &bin : bins()) {
        if (chrom_it == last_chrom || start_it == last_start || end_it == last_end) {
          throw std::runtime_error(
              fmt::format(FMT_STRING("Expected {} bins, found {}"), bins().size(), i));
        }

        if (chromosomes().at(*chrom_it).name() != bin.chrom().name() || *start_it != bin.start() ||
            *end_it != bin.end()) {
          throw std::runtime_error(
              fmt::format(FMT_STRING("Bin #{}: expected {}:{}-{}, found {:ucsc}"), i,
                          chromosomes().at(*chrom_it).name(), *start_it, *end_it, bin));
        }
        ++chrom_it;
        ++start_it;
        ++end_it;
        ++i;
      }
    }

  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Bin table at URI {}/{} is invalid or corrupted: {}"), uri(),
                    group("bins")().getPath(), e.what()));
  }
}

template <typename PixelIt>
inline void File::validate_pixels_before_append(const PixelIt &first_pixel,
                                                const PixelIt &last_pixel) const {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);
  try {
    validate_pixel_type<T>();

    PixelT previous_pixel{};

    std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
      if (pixel.count == T{0}) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) found a pixel of value 0"), pixel));
      }

      if (!chromosomes().contains(pixel.coords.bin1.chrom().id())) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) invalid chromosome id {}"), pixel,
                                             pixel.coords.bin1.chrom().id()));
      }

      if (pixel.coords.bin1.chrom().id() != pixel.coords.bin2.chrom().id() &&
          !chromosomes().contains(pixel.coords.bin2.chrom().id())) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) invalid chromosome id {}"), pixel,
                                             pixel.coords.bin2.chrom().id()));
      }

      if (const auto bin_id = pixel.coords.bin1.id(); bin_id > bins().size()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"),
                        pixel, bin_id));
      }

      if (const auto bin_id = pixel.coords.bin2.id(); bin_id > bins().size()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"),
                        pixel, bin_id));
      }

      if (pixel.coords.bin1.id() > pixel.coords.bin2.id()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}) bin1_id is greater than bin2_id: {} > {}"), pixel,
                        pixel.coords.bin1.id(), pixel.coords.bin2.id()));
      }

      if (!!previous_pixel && previous_pixel.coords >= pixel.coords) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}; {}) pixels are not sorted in ascending order"),
                        previous_pixel.coords, pixel.coords));
      }
      previous_pixel = pixel;
    });

    if (!dataset("pixels/bin1_id").empty()) {
      const auto last_bin1 = dataset("pixels/bin1_id").read_last<std::size_t>();
      const auto last_bin2 = dataset("pixels/bin2_id").read_last<std::size_t>();

      const auto new_bin1 = first_pixel->coords.bin1;
      const auto new_bin2 = first_pixel->coords.bin2;

      if (last_bin1 == new_bin1.id()) {
        if (last_bin2 >= new_bin2.id()) {
          const auto &coord1 = new_bin2;
          const auto coord2 = bins().at(last_bin2);
          throw std::runtime_error(fmt::format(
              FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
        }
      } else if (last_bin1 >= new_bin1.id()) {
        const auto &coord1 = new_bin1;
        const auto coord2 = bins().at(last_bin1);
        throw std::runtime_error(fmt::format(
            FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("pixel validation failed: {}"), e.what()));
  }
}

template <typename PixelIt>
inline void File::validate_thin_pixels_before_append(const PixelIt &first_pixel,
                                                     const PixelIt &last_pixel) const {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);
  validate_pixel_type<T>();

  PixelT previous_pixel{};

  auto pixel_lt_op = [&](const auto &p1, const auto &p2) {
    if (p1.bin1_id != p2.bin1_id) {
      return p1.bin1_id < p2.bin1_id;
    }
    return p1.bin2_id < p2.bin2_id;
  };

  try {
    std::for_each(first_pixel, last_pixel, [&](const ThinPixel<T> &pixel) {
      if (pixel.count == T{0}) {
        throw std::runtime_error("found a pixel of value 0");
      }

      if (pixel.bin1_id >= bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), pixel.bin1_id));
      }

      if (pixel.bin2_id >= bins().size()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("invalid bin id {}: bin maps outside of the bin table"), pixel.bin2_id));
      }

      if (pixel.bin1_id > pixel.bin2_id) {
        throw std::runtime_error(fmt::format(FMT_STRING("bin1_id is greater than bin2_id: {} > {}"),
                                             pixel.bin1_id, pixel.bin2_id));
      }

      if (!!previous_pixel && !pixel_lt_op(previous_pixel, pixel)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}; {}) pixels are not sorted in ascending order"),
                        previous_pixel, pixel));
      }
      previous_pixel = pixel;
    });

    if (!dataset("pixels/bin1_id").empty()) {
      const auto last_bin1 = dataset("pixels/bin1_id").read_last<std::size_t>();
      const auto last_bin2 = dataset("pixels/bin2_id").read_last<std::size_t>();

      const auto new_bin1 = first_pixel->bin1_id;
      const auto new_bin2 = first_pixel->bin2_id;

      if (last_bin1 == new_bin1) {
        if (last_bin2 >= new_bin2) {
          const auto coord1 = bins().at(new_bin2);
          const auto coord2 = bins().at(last_bin2);
          throw std::runtime_error(fmt::format(
              FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
        }
      } else if (last_bin1 >= new_bin1) {
        const auto coord1 = bins().at(new_bin1);
        const auto coord2 = bins().at(last_bin1);
        throw std::runtime_error(fmt::format(
            FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("pixel validation failed: {}"), e.what()));
  }
}

template <typename T1, typename T2>
void assert_holds_alternative(const T2 &buff) {
  if (buff.has_value()) {
    assert(std::holds_alternative<T1>(*buff));
  }
}

template <typename PixelT>
inline void File::validate_pixel_type() const noexcept {
  static_assert(std::is_arithmetic_v<PixelT>);

  if constexpr (std::is_floating_point_v<PixelT>) {
    assert(has_float_pixels());
    assert_holds_alternative<double>(_attrs.sum);
    assert_holds_alternative<double>(_attrs.cis);
  } else {
    assert(has_integral_pixels());
    assert_holds_alternative<std::int64_t>(_attrs.sum);
    assert_holds_alternative<std::int64_t>(_attrs.cis);
  }
}

}  // namespace hictk::cooler
