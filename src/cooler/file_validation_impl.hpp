// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5Exception.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/cooler/validation.hpp"
#include "hictk/fmt.hpp"

namespace hictk::cooler {

inline void File::validate_bins() const {
  try {
    assert(this->_attrs.bin_type == "fixed");
    auto nchroms = this->dataset("bins/chrom").size();
    auto nstarts = this->dataset("bins/start").size();
    auto nends = this->dataset("bins/end").size();
    if (nchroms != nstarts || nchroms != nends) {
      throw std::runtime_error(fmt::format(FMT_STRING("Datasets have inconsistent sizes:\n"
                                                      " - \"bins/chrom\": {}\n"
                                                      " - \"bins/start\": {}\n"
                                                      " - \"bins/end\": {}\n"
                                                      "Expected {}"),
                                           nchroms, nstarts, nends, this->bins().size()));
    }

    const auto &nbins = nchroms;
    if (nbins != this->bins().size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), nchroms));
    }

    auto chrom_it = this->dataset("bins/chrom").begin<std::uint32_t>();
    auto start_it = this->dataset("bins/start").begin<std::uint32_t>();
    auto end_it = this->dataset("bins/end").begin<std::uint32_t>();

    auto last_chrom = this->dataset("bins/chrom").end<std::uint32_t>();
    auto last_start = this->dataset("bins/start").end<std::uint32_t>();
    auto last_end = this->dataset("bins/end").end<std::uint32_t>();

    std::size_t i = 0;
    for (const Bin &bin : this->bins()) {
      if (chrom_it == last_chrom || start_it == last_start || end_it == last_end) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Expected {} bins, found {}"), this->bins().size(), i));
      }

      if (this->chromosomes().at(*chrom_it).name() != bin.chrom().name() ||
          *start_it != bin.start() || *end_it != bin.end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("GenomicInterval #{}: expected {}:{}-{}, found {:ucsc}"), i,
                        this->chromosomes().at(*chrom_it).name(), *start_it, *end_it, bin));
      }
      ++chrom_it;
      ++start_it;
      ++end_it;
      ++i;
    }

  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("GenomicInterval table at URI {}/{} is invalid or corrupted: {}"),
                    this->uri(), this->group("bins")().getPath(), e.what()));
  }
}

template <typename PixelIt>
inline void File::validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);
  try {
    this->validate_pixel_type<T>();

    PixelT previous_pixel{};

    std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
      if (pixel.count == T{0}) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) found a pixel of value 0"), pixel));
      }

      if (!this->chromosomes().contains(pixel.coords.bin1.chrom().id())) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) invalid chromosome id {}"), pixel,
                                             pixel.coords.bin1.chrom().id()));
      }

      if (pixel.coords.bin1.chrom().id() != pixel.coords.bin2.chrom().id() &&
          !this->chromosomes().contains(pixel.coords.bin2.chrom().id())) {
        throw std::runtime_error(fmt::format(FMT_STRING("({}) invalid chromosome id {}"), pixel,
                                             pixel.coords.bin2.chrom().id()));
      }

      if (const auto bin_id = pixel.coords.bin1.id(); bin_id > this->bins().size()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"),
                        pixel, bin_id));
      }

      if (const auto bin_id = pixel.coords.bin2.id(); bin_id > this->bins().size()) {
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

    if (!this->dataset("pixels/bin1_id").empty()) {
      const auto last_bin1 = this->dataset("pixels/bin1_id").read_last<std::size_t>();
      const auto last_bin2 = this->dataset("pixels/bin2_id").read_last<std::size_t>();

      const auto new_bin1 = first_pixel->coords.bin1;
      const auto new_bin2 = first_pixel->coords.bin2;

      if (last_bin1 == new_bin1.id()) {
        if (last_bin2 >= new_bin2.id()) {
          const auto &coord1 = new_bin2;
          const auto coord2 = this->bins().at(last_bin2);
          throw std::runtime_error(fmt::format(
              FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
        }
      } else if (last_bin1 >= new_bin1.id()) {
        const auto &coord1 = new_bin1;
        const auto coord2 = this->bins().at(last_bin1);
        throw std::runtime_error(fmt::format(
            FMT_STRING("new pixel {} is located upstream of pixel {}"), coord1, coord2));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("pixel validation failed: {}"), e.what()));
  }
}

template <typename PixelT>
inline void File::validate_pixel_type() const noexcept {
  static_assert(std::is_arithmetic_v<PixelT>);

  auto assert_holds_alternative = [](const auto &buff, [[maybe_unused]] auto alt) {
    using T [[maybe_unused]] = decltype(alt);
    if (buff.has_value()) {
      assert(std::holds_alternative<T>(*buff));
    }
  };

  if constexpr (std::is_floating_point_v<PixelT>) {
    assert(this->has_float_pixels());
    assert_holds_alternative(this->_attrs.sum, double{});
    assert_holds_alternative(this->_attrs.cis, double{});
  } else {
    assert(this->has_integral_pixels());
    assert_holds_alternative(this->_attrs.sum, std::int64_t{});
    assert_holds_alternative(this->_attrs.cis, std::int64_t{});
  }
}

}  // namespace hictk::cooler
