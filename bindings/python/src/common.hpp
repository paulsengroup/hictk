// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <string>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/pixel.hpp"

namespace hictkpy {

namespace py = pybind11;
using namespace pybind11::literals;

template <typename T>
struct Dynamic1DA {
 private:
  py::array_t<T> _buff{};
  std::int64_t _size{};

 public:
  inline explicit Dynamic1DA(std::size_t size_ = 1000)
      : _buff({static_cast<std::int64_t>(size_)}) {}
  inline void append(T x) {
    if (_buff.size() == _size) {
      grow();
    }
    _buff.mutable_at(_size++) = x;
  }
  inline void grow() { _buff.resize({_buff.size() * 2}); }
  inline void shrink_to_fit() { _buff.resize({_size}); }
  [[nodiscard]] py::array_t<T>&& operator()() noexcept { return std::move(_buff); }
};

template <typename File>
inline py::dict get_chromosomes_from_file(const File& f) {
  py::dict py_chroms{};  // NOLINT
  for (const auto& chrom : f.chromosomes()) {
    const std::string name{chrom.name()};
    py_chroms[name.c_str()] = chrom.size();
  }

  return py_chroms;
}

template <typename File>
inline py::object get_bins_from_file(const File& f) {
  auto pd = py::module::import("pandas");

  std::vector<py::str> chrom_names{};
  Dynamic1DA<std::uint32_t> starts{};
  Dynamic1DA<std::uint32_t> ends{};
  for (const auto& bin : f.bins()) {
    chrom_names.emplace_back(std::string{bin.chrom().name()});
    starts.append(bin.start());
    ends.append(bin.end());
  }

  chrom_names.shrink_to_fit();
  starts.shrink_to_fit();
  ends.shrink_to_fit();

  py::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] = std::move(chrom_names);
  py_bins_dict["start"] = starts();
  py_bins_dict["end"] = ends();

  auto df = pd.attr("DataFrame")(py_bins_dict);
  df.attr("columns") = std::vector<std::string>{"chrom", "start", "end"};
  return df;
}

template <typename PixelIt>
inline py::object pixel_iterators_to_coo(PixelIt first_pixel, PixelIt last_pixel) {
  using N = decltype(first_pixel->count);

  auto pd = py::module::import("pandas");

  Dynamic1DA<std::int64_t> bin1_ids{};
  Dynamic1DA<std::int64_t> bin2_ids{};
  Dynamic1DA<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N>& tp) {
    bin1_ids.append(static_cast<std::int64_t>(tp.bin1_id));
    bin2_ids.append(static_cast<std::int64_t>(tp.bin2_id));
    counts.append(tp.count);
  });

  bin1_ids.shrink_to_fit();
  bin2_ids.shrink_to_fit();
  counts.shrink_to_fit();

  py::dict py_pixels_dict{};  // NOLINT

  py_pixels_dict["bin1_id"] = bin1_ids();
  py_pixels_dict["bin2_id"] = bin2_ids();
  py_pixels_dict["count"] = counts();

  return pd.attr("DataFrame").attr("from_dict")(py_pixels_dict);
}

template <typename PixelIt>
inline py::object pixel_iterators_to_bg2(const hictk::BinTable& bins, PixelIt first_pixel,
                                         PixelIt last_pixel) {
  using N = decltype(first_pixel->count);

  auto pd = py::module::import("pandas");

  std::vector<py::str> chrom_names1{};
  Dynamic1DA<std::int32_t> starts1{};
  Dynamic1DA<std::int32_t> ends1{};
  std::vector<py::str> chrom_names2{};
  Dynamic1DA<std::int32_t> starts2{};
  Dynamic1DA<std::int32_t> ends2{};
  Dynamic1DA<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const hictk::ThinPixel<N>& tp) {
    const hictk::Pixel<N> p{bins, tp};

    chrom_names1.emplace_back(p.coords.bin1.chrom().name());
    starts1.append(static_cast<std::int32_t>(p.coords.bin1.start()));
    ends1.append(static_cast<std::int32_t>(p.coords.bin1.end()));

    chrom_names2.emplace_back(p.coords.bin2.chrom().name());
    starts2.append(static_cast<std::int32_t>(p.coords.bin2.start()));
    ends2.append(static_cast<std::int32_t>(p.coords.bin2.end()));

    counts.append(p.count);
  });

  starts1.shrink_to_fit();
  ends1.shrink_to_fit();
  starts2.shrink_to_fit();
  ends2.shrink_to_fit();
  counts.shrink_to_fit();

  py::dict py_pixels_dict{};  // NOLINT

  py_pixels_dict["chrom1"] = py::array(py::cast(chrom_names1));
  py_pixels_dict["start1"] = starts1();
  py_pixels_dict["end1"] = ends1();
  py_pixels_dict["chrom2"] = py::array(py::cast(chrom_names2));
  py_pixels_dict["start2"] = starts2();
  py_pixels_dict["end2"] = ends2();

  py_pixels_dict["count"] = counts();

  return pd.attr("DataFrame").attr("from_dict")(py_pixels_dict);
}

template <typename PixelIt>
static py::object pixel_iterators_to_df(const hictk::BinTable& bins, PixelIt first_pixel,
                                        PixelIt last_pixel, bool join) {
  if (join) {
    return pixel_iterators_to_bg2(bins, first_pixel, last_pixel);
  }
  return pixel_iterators_to_coo(first_pixel, last_pixel);
}

template <typename File>
inline py::object file_fetch_all(File& f, std::string_view normalization,
                                 std::string_view count_type, bool join) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }

  if (normalization != "NONE") {
    count_type = "float";
  }

  auto sel = f.fetch(hictk::balancing::Method{normalization});
  if (count_type == "int") {
    return pixel_iterators_to_df(f.bins(), sel.template begin<std::int32_t>(),
                                 sel.template end<std::int32_t>(), join);
  }
  return pixel_iterators_to_df(f.bins(), sel.template begin<double>(), sel.template end<double>(),
                               join);
}

template <typename File>
inline py::object file_fetch(const File& f, std::string_view range1, std::string_view range2,
                             std::string_view normalization, std::string_view count_type, bool join,
                             std::string_view query_type) {
  if (range1.empty()) {
    return file_fetch_all(f, normalization, count_type, join);
  }
  if (normalization != "NONE") {
    count_type = "float";
  }

  const auto qt =
      query_type == "UCSC" ? hictk::GenomicInterval::Type::UCSC : hictk::GenomicInterval::Type::BED;

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hictk::balancing::Method(normalization), qt)
                 : f.fetch(range1, range2, hictk::balancing::Method(normalization), qt);

  if (count_type == "int") {
    return pixel_iterators_to_df(f.bins(), sel.template begin<std::int32_t>(),
                                 sel.template end<std::int32_t>(), join);
  }
  return pixel_iterators_to_df(f.bins(), sel.template begin<double>(), sel.template end<double>(),
                               join);
}

}  // namespace hictkpy
