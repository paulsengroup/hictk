// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "./validate.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/file_attributes_formatting.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

static void update_status_table(const cooler::utils::ValidationStatusCooler& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("is_valid_cooler", status.is_cooler);
}

[[nodiscard]] static bool validate_bin_table_shape(const cooler::File& clr) {
  const auto& chroms = clr.dataset("bins/chrom");
  const auto& starts = clr.dataset("bins/start");
  const auto& ends = clr.dataset("bins/end");

  const auto expected_num_bins = clr.bins().size();
  return chroms.size() == expected_num_bins && starts.size() == expected_num_bins &&
         ends.size() == expected_num_bins;
}

[[nodiscard]] static bool validate_bins_dtypes(const cooler::File& clr) {
  try {
    std::ignore = *clr.dataset("bins/chrom").begin<std::string>();
    std::ignore = *clr.dataset("bins/start").begin<std::int32_t>();
    std::ignore = *clr.dataset("bins/end").begin<std::int32_t>();
    return true;
  } catch (...) {
    return false;
  }
}

[[nodiscard]] static std::size_t count_invalid_bins(const cooler::File& clr) {
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

  return num_invalid_bins;
}

static bool check_bin_table(const cooler::File& clr, toml::table& status) {
  const auto size_ok = validate_bin_table_shape(clr);
  status.insert("bin_table_shape_ok", size_ok);
  if (!size_ok) {
    return false;
  }

  const auto dtypes_ok = validate_bins_dtypes(clr);
  status.insert("bin_table_dtypes_ok", dtypes_ok);
  if (!dtypes_ok) {
    return false;
  }

  const auto num_invalid_bins = count_invalid_bins(clr);
  status.insert("bin_table_num_invalid_bins", static_cast<std::int64_t>(num_invalid_bins));
  return num_invalid_bins == 0;
}

template <typename N>
static void check_duplicate_pixel([[maybe_unused]] const cooler::File& clr,
                                  const ThinPixel<N>& prev_pixel, const ThinPixel<N>& pixel,
                                  std::size_t i, std::optional<std::string>& status) {
  if (HICTK_UNLIKELY(status.has_value())) {
    return;
  }

  if (HICTK_UNLIKELY(pixel.bin1_id == prev_pixel.bin1_id && pixel.bin2_id == prev_pixel.bin2_id)) {
    status = fmt::format(
        FMT_STRING("pixel #{} and #{} have the same coordinates (bin1_id={} and bin2_id={})"),
        i - 1, i, pixel.bin1_id, pixel.bin2_id);
    SPDLOG_DEBUG(FMT_STRING("{}: {}"), clr.uri(), *status);
  }
}

template <typename N>
static void check_matrix_symmetry([[maybe_unused]] const cooler::File& clr,
                                  const ThinPixel<N>& pixel, bool file_is_symmetric_upper,
                                  std::size_t i, std::optional<std::string>& status) {
  if (HICTK_UNLIKELY(!file_is_symmetric_upper || status.has_value())) {
    return;
  }

  if (pixel.bin1_id > pixel.bin2_id) {
    status = fmt::format(
        FMT_STRING("pixel #{} (bin1_id={} bin2_id={}) overlaps with the lower-triangular matrix"),
        i, pixel.bin1_id, pixel.bin2_id);
    SPDLOG_DEBUG(FMT_STRING("{}: {}"), clr.uri(), *status);
  }
}

template <typename N>
static void check_pixels_are_sorted([[maybe_unused]] const cooler::File& clr,
                                    const ThinPixel<N>& prev_pixel, const ThinPixel<N>& pixel,
                                    std::size_t i, std::optional<std::string>& status) {
  if (HICTK_UNLIKELY(status.has_value())) {
    return;
  }

  if (HICTK_UNLIKELY(prev_pixel > pixel)) {
    status = fmt::format(
        FMT_STRING("pixel #{} and #{} are not sorted in ascending order: {}:{} > {}:{}"), i - 1, i,
        prev_pixel.bin1_id, prev_pixel.bin2_id, pixel.bin1_id, pixel.bin2_id);
    SPDLOG_DEBUG(FMT_STRING("{}: {}"), clr.uri(), *status);
  }
}

template <typename N>
static void check_pixel_count([[maybe_unused]] const cooler::File& clr, const ThinPixel<N>& pixel,
                              std::size_t i, std::optional<std::string>& status) {
  if (HICTK_UNLIKELY(status.has_value())) {
    return;
  }
  if (HICTK_UNLIKELY(pixel.count == 0)) {
    status = fmt::format(FMT_STRING("pixel #{} has an invalid count {}:{}={}"), i, pixel.bin1_id,
                         pixel.bin2_id, pixel.count);
    SPDLOG_DEBUG(FMT_STRING("{}: {}"), clr.uri(), *status);
  }
}

template <typename N>
static bool check_pixels(const cooler::File& clr, toml::table& status) {
  SPDLOG_DEBUG(FMT_STRING("{}: validating pixels..."), clr.uri());
  auto first = clr.begin<N>();
  const auto last = clr.end<N>();

  if (first == last) {
    return true;
  }

  const auto symmetric_upper =
      clr.attributes().storage_mode.value_or("symmetric-upper") == "symmetric-upper";

  auto prev_pixel = *first++;
  std::optional<std::string> dupl_pixel_status{};
  std::optional<std::string> sorted_status{};
  std::optional<std::string> count_status{};
  std::optional<std::string> symmetry_status{};

  auto early_return = [&]() {
    return symmetry_status.has_value() && dupl_pixel_status.has_value() &&
           count_status.has_value() && sorted_status.has_value();
  };

  check_pixel_count(clr, prev_pixel, 1, count_status);
  check_matrix_symmetry(clr, prev_pixel, symmetric_upper, 1, symmetry_status);

  for (std::size_t i = 2; first != last; ++i) {
    auto pixel = *first;
    check_pixel_count(clr, pixel, i, count_status);
    check_matrix_symmetry(clr, pixel, symmetric_upper, i, symmetry_status);
    check_pixels_are_sorted(clr, prev_pixel, pixel, i, sorted_status);
    check_duplicate_pixel(clr, prev_pixel, pixel, i, dupl_pixel_status);
    if (HICTK_UNLIKELY(early_return())) {
      break;
    }
    std::swap(prev_pixel, pixel);
    std::ignore = ++first;
  }

  if (sorted_status.has_value()) {
    status.insert("pixels_are_sorted", *sorted_status);
  } else {
    status.insert("pixels_are_sorted", true);
  }
  if (!symmetric_upper) {
    status.insert("pixels_are_symmetric_upper", "not_checked");
  } else if (symmetry_status.has_value()) {
    status.insert("pixels_are_symmetric_upper", *symmetry_status);
  } else {
    status.insert("pixels_are_symmetric_upper", true);
  }
  if (dupl_pixel_status.has_value()) {
    status.insert("pixels_are_unique", *dupl_pixel_status);
  } else {
    status.insert("pixels_are_unique", true);
  }
  if (count_status.has_value()) {
    status.insert("pixels_have_valid_counts", *count_status);
  } else {
    status.insert("pixels_have_valid_counts", true);
  }

  return !sorted_status.has_value() && !symmetry_status.has_value() &&
         !dupl_pixel_status.has_value() && !count_status.has_value();
}

static bool check_pixels(const cooler::File& clr, toml::table& status) {
  if (clr.has_float_pixels()) {
    return check_pixels<double>(clr, status);
  }
  if (clr.has_unsigned_pixels()) {
    return check_pixels<std::uint64_t>(clr, status);
  }
  return check_pixels<std::int64_t>(clr, status);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
std::pair<int, toml::table> validate_cooler(std::string_view path, bool validate_index,
                                            bool validate_pixels) {
  int return_code = 0;
  toml::table status;

  update_status_table(cooler::utils::is_cooler(path), status);
  auto is_cooler = **status.get_as<bool>("is_valid_cooler");

  std::optional<cooler::File> clr{};
  if (is_cooler) {
    try {
      clr = cooler::File(path);
    } catch (...) {
      is_cooler = false;
      status.insert_or_assign("is_valid_cooler", is_cooler);
      return_code = 1;
    }
  }

  if (!!clr) {
    assert(clr.has_value());
    check_bin_table(*clr, status);
  } else {
    status.insert("bin_table_shape_ok", "not_checked");
    status.insert("bin_table_dtypes_ok", "not_checked");
    status.insert("bin_table_num_invalid_bins", "not_checked");
  }

  if (!!clr && validate_index) {
    std::string buff{};
    try {
      const auto index_ok = cooler::utils::index_is_valid(path, buff);
      if (index_ok) {
        status.insert("index_is_valid", true);
      } else {
        assert(!buff.empty());
        return_code = 1;
        status.insert("index_is_valid", buff);
      }
    } catch (const std::exception& e) {
      if (const auto storage_mode_ok =
              !internal::starts_with(e.what(),
                                     "validating the index of Coolers with storage-mode") &&
              !internal::ends_with(e.what(), "is not supported");
          storage_mode_ok) {
        throw;
      }
      SPDLOG_WARN(FMT_STRING("{}"), e.what());
      status.insert("index_is_valid", "not_checked");
    }
  } else {
    status.insert("index_is_valid", "not_checked");
  }

  if (!!clr && validate_pixels) {
    if (!check_pixels(*clr, status)) {
      return_code = 1;
    }
  } else {
    status.insert("pixels_are_sorted", "not_checked");
    status.insert("pixels_are_symmetric_upper", "not_checked");
    status.insert("pixels_are_unique", "not_checked");
    status.insert("pixels_have_valid_counts", "not_checked");
  }

  if (return_code != 0) {
    status.insert_or_assign("is_valid_cooler", false);
  }

  return std::make_pair(return_code, status);
}

}  // namespace hictk::tools
