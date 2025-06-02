// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "./common.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/string.hpp"
#include "hictk/tools/compressed_io.hpp"

namespace hictk::tools {

struct Header {
  std::string assembly{"unknown"};
  std::optional<Reference> chromosomes{};
};

class PixelParser {
  io::CompressedReader _reader{};
  Format _format{Format::_4DN};
  std::string _strbuff{};
  BinTable _bins{};
  std::string _assembly{};
  bool _drop_unknown_chroms{false};
  std::size_t _num_dropped_records{};

 public:
  PixelParser() = default;
  PixelParser(const std::filesystem::path& path_, std::uint32_t resolution, Format format_,
              std::string_view assembly = "unknown", bool drop_unknown_chroms = false);
  PixelParser(const std::filesystem::path& path_, BinTable bins_, Format format_,
              std::string_view assembly_ = "unknown", bool drop_unknown_chroms = false);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string_view assembly() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;

  template <typename N>
  [[nodiscard]] ThinPixel<N> next_pixel(std::int64_t offset) {
    ThinPixel<N> p{};
    std::ignore = next_pixel(p, offset);
    return p;
  }

  template <typename N>
  [[nodiscard]] bool next_pixel(ThinPixel<N>& buff, std::int64_t offset) {
    auto handle_chrom_not_found = [&](std::string_view msg) {
      if (!_drop_unknown_chroms) {
        throw;
      }

      const auto pos = msg.rfind("chromosome \"");
      if (pos == std::string_view::npos) {
        throw;
      }

      msg.remove_prefix(pos);
      const auto chrom_not_found = internal::ends_with(msg, "\" not found");
      if (chrom_not_found) {
        ++_num_dropped_records;
        std::ignore = getline();
        return;
      }

      throw;
    };

    while (HICTK_LIKELY(!_strbuff.empty())) {
      try {
        switch (_format) {
          case Format::COO:
            buff = ThinPixel<N>::from_coo(bins(), _strbuff, offset);
            break;
          case Format::BG2:
            buff = Pixel<N>::from_bg2(bins(), _strbuff, offset).to_thin();
            break;
          case Format::VP:
            buff = Pixel<N>::from_validpair(bins(), _strbuff, offset).to_thin();
            break;
          case Format::_4DN:
            buff = Pixel<N>::from_4dn_pairs(bins(), _strbuff, offset).to_thin();
            break;
        }
        std::ignore = getline();
        return true;
      } catch (const std::out_of_range& e) {
        handle_chrom_not_found(e.what());
      } catch (const std::runtime_error& e) {
        handle_chrom_not_found(e.what());
      }
    }

    // EOF
    assert(_strbuff.empty());
    buff.bin1_id = ThinPixel<N>::null_id;
    buff.bin2_id = ThinPixel<N>::null_id;
    buff.count = 0;
    return false;
  }

 private:
  [[nodiscard]] bool getline(char sep = '\n');

  [[nodiscard]] static std::string_view remove_prefix(std::string_view s, std::string_view prefix,
                                                      bool strip_whitespaces = true) noexcept;
  [[nodiscard]] static std::pair<std::string, std::uint32_t> parse_chromsize(std::string_view line);

  [[nodiscard]] Header parse_header(Format format_);
};

[[nodiscard]] PixelParser init_pixel_parser(Format format,
                                            const std::filesystem::path& path_to_interactions,
                                            const std::filesystem::path& path_to_chrom_sizes,
                                            const std::filesystem::path& path_to_bins,
                                            std::uint32_t resolution, std::string_view assembly,
                                            bool drop_unknown_chroms);

}  // namespace hictk::tools
