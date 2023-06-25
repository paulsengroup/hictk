// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <string_view>

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

template <typename N>
[[nodiscard]] static Pixel<N> parse_pixel(const BinTable& bins, std::string_view line) {
  auto next_token = [&]() {
    assert(!line.empty());
    const auto pos = line.find('\t');
    auto tok = line.substr(0, pos);
    line.remove_prefix(pos + 1);
    return tok;
  };

  const auto chrom1 = next_token();
  const auto start1 = internal::parse_numeric_or_throw<std::uint32_t>(next_token());
  std::ignore = next_token();

  const auto chrom2 = next_token();
  const auto start2 = internal::parse_numeric_or_throw<std::uint32_t>(next_token());
  std::ignore = next_token();

  const auto count = internal::parse_numeric_or_throw<N>(next_token());
  return Pixel<N>{{bins.at(chrom1, start1), bins.at(chrom2, start2)}, count};
}

template <typename N>
static void read_batch(const BinTable& bins, std::vector<Pixel<N>>& buffer) {
  buffer.clear();
  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (buffer.size() == buffer.capacity()) {
        return;
      }
      buffer.emplace_back(parse_pixel<N>(bins, line));
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }
}

template <typename N>
static std::string ingest_pixels_sorted(hictk::cooler::File&& clr,
                                        std::size_t batch_size = 500'000) {
  const auto& bins = clr.bins();

  std::vector<Pixel<N>> buffer(batch_size);

  buffer.clear();
  std::string line;
  std::size_t i = 0;
  try {
    for (; std::getline(std::cin, line); ++i) {
      if (buffer.size() == buffer.capacity()) {
        clr.append_pixels(buffer.begin(), buffer.end(), true);
        buffer.clear();
      }
      buffer.emplace_back(parse_pixel<N>(bins, line));
    }
    clr.append_pixels(buffer.begin(), buffer.end(), true);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while processing line(s) {}-{}: {}"),
                    i - buffer.size(), i, e.what()));
  }
  return clr.uri();
}

template <typename N>
static std::string ingest_pixels_unsorted(hictk::cooler::File&& clr,
                                          std::vector<Pixel<N>>& write_buffer) {
  assert(write_buffer.capacity() != 0);

  if (std::cin.eof()) {
    return "";
  }

  read_batch(clr.bins(), write_buffer);
  if (write_buffer.empty()) {
    return "";
  }

  std::sort(write_buffer.begin(), write_buffer.end());
  clr.append_pixels(write_buffer.begin(), write_buffer.end());

  return clr.uri();
}

void load_subcmd(const LoadConfig& c) {
  auto chroms = Reference(c.path_to_chrom_sizes);
  if (c.assume_sorted) {
    c.count_as_float
        ? ingest_pixels_sorted<double>(
              hictk::cooler::File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force))
        : ingest_pixels_sorted<std::int32_t>(hictk::cooler::File::create_new_cooler<std::int32_t>(
              c.uri, chroms, c.bin_size, c.force));
    return;
  }

  std::vector<std::string> uris{};
  std::vector<Pixel<std::int32_t>> write_buffer_int{};
  std::vector<Pixel<double>> write_buffer_fp{};

  const internal::TmpDir tmpdir{};
  for (std::size_t i = 0; true; ++i) {
    const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
    if (c.count_as_float) {
      write_buffer_fp.reserve(c.batch_size);
      uris.emplace_back(ingest_pixels_unsorted(hictk::cooler::File::create_new_cooler<double>(
                                                   tmp_uri.string(), chroms, c.bin_size, c.force),
                                               write_buffer_fp));
    } else {
      write_buffer_int.reserve(c.batch_size);
      uris.emplace_back(ingest_pixels_unsorted(hictk::cooler::File::create_new_cooler<std::int32_t>(
                                                   tmp_uri.string(), chroms, c.bin_size, c.force),
                                               write_buffer_int));
    }
    if (uris.back().empty()) {
      uris.pop_back();
      break;
    }
    spdlog::info(FMT_STRING("Done writing to tmp file {}..."), tmp_uri);
  }

  if (c.force) {
    std::filesystem::remove(c.uri);
  }

  if (uris.size() == 1) {
    spdlog::info(FMT_STRING("Moving temporary file to {}..."), c.uri);
    std::filesystem::copy_file(uris.front(), c.uri);
    return;
  }

  MergeConfig cm{uris, c.uri};
  cm.force = c.force;
  cm.verbosity = c.verbosity;

  spdlog::info(FMT_STRING("Merging {} intermediate files..."), uris.size());
  hictk::cooler::utils::merge(uris.begin(), uris.end(), c.uri, c.force, 500'000, c.verbosity != 3);
}

}  // namespace hictk::tools
