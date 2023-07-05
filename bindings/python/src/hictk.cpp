// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"

namespace py = pybind11;
using namespace hictk;

static bool is_cooler(std::string_view uri) { return bool(cooler::utils::is_cooler(uri)); }

static cooler::File cooler_ctor(std::string_view uri) { return cooler::File::open_read_only(uri); }

static cooler::File cooler_ctor(std::string_view uri, const py::dict& py_chroms,
                                std::uint32_t bin_size, bool overwrite_if_exists = false) {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  for (auto it : py_chroms) {
    chrom_names.push_back(py::cast<std::string>(it.first));
    chrom_sizes.push_back(py::cast<std::uint32_t>(it.second));
  }
  Reference chroms(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
  return cooler::File::create_new_cooler(uri, chroms, bin_size, overwrite_if_exists);
}

static py::dict cooler_chromosomes(const cooler::File& clr) {
  py::dict py_chroms{};
  for (const auto& chrom : clr.chromosomes()) {
    const std::string name{chrom.name()};
    py_chroms[name.c_str()] = chrom.size();
  }

  return py_chroms;
}

static py::object cooler_bins(const cooler::File& clr) {
  auto pd = py::module::import("pandas");

  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> starts{};
  std::vector<std::uint32_t> ends{};
  for (const auto& bin : clr.bins()) {
    chrom_names.emplace_back(bin.chrom().name());
    starts.emplace_back(bin.start());
    ends.emplace_back(bin.end());
  }

  py::dict py_bins_dict{};

  py_bins_dict["chrom"] = std::move(chrom_names);
  py_bins_dict["start"] = std::move(starts);
  py_bins_dict["end"] = std::move(ends);

  auto df = pd.attr("DataFrame")(py_bins_dict);
  df.attr("columns") = std::vector<std::string>{"chrom", "start", "end"};
  return df;
}

template <typename PixelIt>
static py::object pixel_iterators_to_df(PixelIt first_pixel, PixelIt last_pixel) {
  using N = decltype(first_pixel->count);

  auto pd = py::module::import("pandas");

  std::vector<std::string> chrom_names1{};
  std::vector<std::int64_t> starts1{};
  std::vector<std::int64_t> ends1{};

  std::vector<std::string> chrom_names2{};
  std::vector<std::int64_t> starts2{};
  std::vector<std::int64_t> ends2{};
  std::vector<N> counts{};

  std::for_each(first_pixel, last_pixel, [&](const Pixel<N>& p) {
    chrom_names1.emplace_back(p.coords.bin1.chrom().name());
    starts1.emplace_back(p.coords.bin1.start());
    ends1.emplace_back(p.coords.bin1.end());

    chrom_names2.emplace_back(p.coords.bin2.chrom().name());
    starts2.emplace_back(p.coords.bin2.start());
    ends2.emplace_back(p.coords.bin2.end());

    counts.push_back(p.count);
  });

  py::dict py_pixels_dict{};

  py_pixels_dict["chrom1"] = std::move(chrom_names1);
  py_pixels_dict["start1"] = std::move(starts1);
  py_pixels_dict["end1"] = std::move(ends1);

  py_pixels_dict["chrom2"] = std::move(chrom_names2);
  py_pixels_dict["start2"] = std::move(starts2);
  py_pixels_dict["end2"] = std::move(ends2);

  py_pixels_dict["count"] = std::move(counts);

  auto df = pd.attr("DataFrame")(py_pixels_dict);
  df.attr("columns") =
      std::vector<std::string>{"chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"};
  return df;
}

static py::object cooler_fetch(const cooler::File& clr, std::string_view range1,
                               std::string_view range2, std::string_view count_type,
                               std::string_view query_type) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }
  const auto qt =
      query_type == "UCSC" ? cooler::File::QUERY_TYPE::UCSC : cooler::File::QUERY_TYPE::BED;

  if (count_type == "int") {
    auto sel = range2.empty() || range1 == range2 ? clr.fetch<std::int32_t>(range1, qt)
                                                  : clr.fetch<std::int32_t>(range1, range2, qt);
    return pixel_iterators_to_df(sel.begin(), sel.end());
  }
  auto sel = range2.empty() || range1 == range2 ? clr.fetch<double>(range1, qt)
                                                : clr.fetch<double>(range1, range2, qt);
  return pixel_iterators_to_df(sel.begin(), sel.end());
}

namespace internal {
template <typename N>
static py::array_t<N> cooler_fetch_matrix(const cooler::File& clr, std::string_view range1,
                                          std::string_view range2, std::string_view query_type) {
  auto np = py::module::import("numpy");

  const auto qt =
      query_type == "UCSC" ? cooler::File::QUERY_TYPE::UCSC : cooler::File::QUERY_TYPE::BED;
  auto sel = range2.empty() || range1 == range2 ? clr.fetch<N>(range1, qt)
                                                : clr.fetch<N>(range1, range2, qt);

  const auto bin11 = sel.coord1().bin1.id();
  const auto bin12 = sel.coord1().bin2.id() + 1;

  const auto bin21 = sel.coord2().bin1.id();
  const auto bin22 = sel.coord2().bin2.id() + 1;

  const auto rows = bin12 - bin11;
  const auto cols = bin22 - bin21;

  py::array_t<N> m(static_cast<ssize_t>(rows * cols));
  auto* m_ptr = m.mutable_data();

  for (const auto& p : sel) {
    const auto i1 = p.coords.bin1.id() - bin11;
    const auto i2 = p.coords.bin2.id() - bin21;

    const auto i = (i1 * cols) + i2;
    const auto j = (i2 * cols) + i1;

    m_ptr[i] = p.count;
    m_ptr[j] = p.count;
  }

  return m.reshape({rows, cols});
}
}  // namespace internal

static py::object cooler_fetch_matrix(const cooler::File& clr, std::string_view range1,
                                      std::string_view range2, std::string_view count_type,
                                      std::string_view query_type) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }
  if (count_type == "int") {
    return ::internal::cooler_fetch_matrix<std::int32_t>(clr, range1, range2, query_type);
  }
  return ::internal::cooler_fetch_matrix<double>(clr, range1, range2, query_type);
}

PYBIND11_MODULE(hictkpy, m) {
  [[maybe_unused]] auto np = py::module::import("numpy");
  [[maybe_unused]] auto pd = py::module::import("pandas");

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files";
  auto cooler = m.def_submodule("cooler");
  auto cooler_utils = cooler.def_submodule("utils");

  cooler_utils.def("is_cooler", &is_cooler, "test whether path points to a cooler file");

  auto file =
      py::class_<cooler::File>(cooler, "File")
          .def(py::init(py::overload_cast<std::string_view>(cooler_ctor)), py::arg("uri"))
          .def(py::init(py::overload_cast<std::string_view, const py::dict&, std::uint32_t, bool>(
                   cooler_ctor)),
               py::arg("uri"), py::arg("chromosomes"), py::arg("bin_size"),
               py::arg("overwrite_if_exists"))
          .def("open_read_only", &cooler::File::open_read_only);

  file.def("uri", &cooler::File::uri);
  file.def("hdf5_path", &cooler::File::hdf5_path);
  file.def("path", &cooler::File::path);

  file.def("bin_size", &cooler::File::bin_size);
  file.def("chromosomes", &cooler_chromosomes);
  file.def("bins", &cooler_bins);

  file.def("fetch", &cooler_fetch, py::arg("range1"), py::arg("range2") = "",
           py::arg("count_type") = "int", py::arg("query_type") = "UCSC");
  file.def("fetch_matrix", &cooler_fetch_matrix, py::arg("range1"), py::arg("range2") = "",
           py::arg("count_type") = "int", py::arg("query_type") = "UCSC");
}
