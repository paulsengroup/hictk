// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"

namespace py = pybind11;
using namespace pybind11::literals;
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
  const Reference chroms(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
  return cooler::File::create_new_cooler(uri, chroms, bin_size, overwrite_if_exists);
}

template <typename File>
static py::dict get_chromosomes_from_file(const File& f) {
  py::dict py_chroms{};  // NOLINT
  for (const auto& chrom : f.chromosomes()) {
    const std::string name{chrom.name()};
    py_chroms[name.c_str()] = chrom.size();
  }

  return py_chroms;
}

template <typename File>
static py::object get_bins_from_file(const File& f) {
  auto pd = py::module::import("pandas");

  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> starts{};
  std::vector<std::uint32_t> ends{};
  for (const auto& bin : f.bins()) {
    chrom_names.emplace_back(bin.chrom().name());
    starts.emplace_back(bin.start());
    ends.emplace_back(bin.end());
  }

  py::dict py_bins_dict{};  // NOLINT

  py_bins_dict["chrom"] = std::move(chrom_names);
  py_bins_dict["start"] = std::move(starts);
  py_bins_dict["end"] = std::move(ends);

  auto df = pd.attr("DataFrame")(py_bins_dict);
  df.attr("columns") = std::vector<std::string>{"chrom", "start", "end"};
  return df;
}

template <typename T>
struct Dynamic1DA {
 private:
  py::array_t<T> _buff{};
  std::int64_t _size{};

 public:
  explicit Dynamic1DA(std::size_t size_ = 1000) : _buff({static_cast<std::int64_t>(size_)}) {}
  void append(T x) {
    if (_buff.size() == _size) {
      grow();
    }
    _buff.mutable_at(_size++) = x;
  }
  void grow() { _buff.resize({_buff.size() * 2}); }
  void shrink_to_fit() { _buff.resize({_size}); }
  [[nodiscard]] py::array_t<T>&& operator()() noexcept { return std::move(_buff); }
};

template <typename PixelIt>
static py::object pixel_iterators_to_df(const BinTable& bins, PixelIt first_pixel,
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

  std::for_each(first_pixel, last_pixel, [&](const ThinPixel<N>& tp) {
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

static py::object cooler_fetch(const cooler::File& clr, std::string_view range1,
                               std::string_view range2, std::string_view normalization,
                               std::string_view count_type, std::string_view query_type) {
  if (count_type != "int" && count_type != "float") {
    throw std::runtime_error("invalid count type. Allowed types: int, float.");
  }

  const auto weights = clr.read_weights(normalization);
  if (weights) {
    count_type = "float";
  }

  const auto qt =
      query_type == "UCSC" ? cooler::File::QUERY_TYPE::UCSC : cooler::File::QUERY_TYPE::BED;

  if (count_type == "int") {
    auto sel = range2.empty() || range1 == range2 ? clr.fetch(range1, weights, qt)
                                                  : clr.fetch(range1, range2, weights, qt);
    return pixel_iterators_to_df(clr.bins(), sel.begin<std::int32_t>(), sel.end<std::int32_t>());
  }
  auto sel = range2.empty() || range1 == range2 ? clr.fetch(range1, weights, qt)
                                                : clr.fetch(range1, range2, weights, qt);
  return pixel_iterators_to_df(clr.bins(), sel.begin<double>(), sel.end<double>());
}

[[nodiscard]] hic::HiCFile hic_ctor(std::string_view path, std::int32_t resolution,
                                    std::string_view matrix_type, std::string_view matrix_unit) {
  return hic::HiCFile{std::string{path}, static_cast<std::uint32_t>(resolution),
                      hic::ParseMatrixTypeStr(std::string{matrix_type}),
                      hic::ParseUnitStr(std::string{matrix_unit})};
}

static py::object hic_fetch(const hic::HiCFile& f, std::string_view range1, std::string_view range2,
                            std::string_view normalization, std::string_view query_type) {
  const auto qt =
      query_type == "UCSC" ? cooler::File::QUERY_TYPE::UCSC : cooler::File::QUERY_TYPE::BED;

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hic::ParseNormStr(std::string{normalization}), qt)
                 : f.fetch(range1, range2, hic::ParseNormStr(std::string{normalization}), qt);
  if (normalization == "NONE") {
    return pixel_iterators_to_df(f.bins(), sel.begin<std::int32_t>(), sel.end<std::int32_t>());
  }
  return pixel_iterators_to_df(f.bins(), sel.begin<double>(), sel.end<double>());
}

PYBIND11_MODULE(hictkpy, m) {
  [[maybe_unused]] auto np = py::module::import("numpy");
  [[maybe_unused]] auto pd = py::module::import("pandas");
  m.attr("__version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files";
  auto cooler = m.def_submodule("cooler");
  auto cooler_utils = cooler.def_submodule("utils");

  cooler_utils.def("is_cooler", &is_cooler, "test whether path points to a cooler file");

  auto cooler_file =
      py::class_<cooler::File>(cooler, "File")
          .def(py::init(py::overload_cast<std::string_view>(cooler_ctor)), py::arg("uri"))
          .def(py::init(py::overload_cast<std::string_view, const py::dict&, std::uint32_t, bool>(
                   cooler_ctor)),
               py::arg("uri"), py::arg("chromosomes"), py::arg("bin_size"),
               py::arg("overwrite_if_exists"))
          .def("open_read_only", &cooler::File::open_read_only);

  cooler_file.def("uri", &cooler::File::uri);
  cooler_file.def("hdf5_path", &cooler::File::hdf5_path);
  cooler_file.def("path", &cooler::File::path);

  cooler_file.def("bin_size", &cooler::File::bin_size);
  cooler_file.def("chromosomes", &get_chromosomes_from_file<cooler::File>);
  cooler_file.def("bins", &get_bins_from_file<cooler::File>);

  cooler_file.def("fetch", &cooler_fetch, py::arg("range1"), py::arg("range2") = "",
                  py::arg("normalization") = "NONE", py::arg("count_type") = "int",
                  py::arg("query_type") = "UCSC");

  auto hic = m.def_submodule("hic");
  auto hic_utils = hic.def_submodule("utils");

  hic_utils.def("is_hic_file", &hic::utils::is_hic_file, "test whether path points to a .hic file");

  auto hic_file = py::class_<hic::HiCFile>(hic, "File")
                      .def(py::init(&hic_ctor), py::arg("path"), py::arg("resolution"),
                           py::arg("matrix_type") = "observed", py::arg("matrix_unit") = "BP");

  hic_file.def("path", &hic::HiCFile::url);
  hic_file.def("name", &hic::HiCFile::name);
  hic_file.def("version", &hic::HiCFile::version);

  hic_file.def("bin_size", &hic::HiCFile::resolution);
  hic_file.def("chromosomes", &get_chromosomes_from_file<hic::HiCFile>);
  hic_file.def("bins", &get_bins_from_file<hic::HiCFile>);

  hic_file.def("fetch", &hic_fetch, py::arg("range1"), py::arg("range2") = "",
               py::arg("normalization") = "NONE", py::arg("query_type") = "UCSC");
}
