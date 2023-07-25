// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <pybind11/pybind11.h>

#include "./common.hpp"
#include "./hictkpy_cooler.hpp"
#include "./hictkpy_file.hpp"
#include "./hictkpy_hic.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"

namespace hictkpy {

static pybind11::module_ declare_cooler_submodule(pybind11::module_& m) {
  auto cooler = m.def_submodule("cooler");
  auto cooler_utils = cooler.def_submodule("utils");

  cooler_utils.def("is_cooler", &cooler::is_cooler, "test whether path points to a cooler file");

  auto cooler_file =
      py::class_<hictk::cooler::File>(cooler, "File")
          .def(py::init(py::overload_cast<std::string_view>(cooler::file_ctor)), py::arg("uri"))
          .def(py::init(py::overload_cast<std::string_view, const py::dict&, std::uint32_t, bool>(
                   cooler::file_ctor)),
               py::arg("uri"), py::arg("chromosomes"), py::arg("bin_size"),
               py::arg("overwrite_if_exists"));

  cooler_file.def("uri", &hictk::cooler::File::uri);
  cooler_file.def("hdf5_path", &hictk::cooler::File::hdf5_path);
  cooler_file.def("path", &hictk::cooler::File::path);

  cooler_file.def("bin_size", &hictk::cooler::File::bin_size);
  cooler_file.def("nbins", &hictk::cooler::File::nbins);
  cooler_file.def("nchroms", &hictk::cooler::File::nchroms);
  cooler_file.def("nnz", &hictk::cooler::File::nnz);

  cooler_file.def("chromosomes", &get_chromosomes_from_file<hictk::cooler::File>);
  cooler_file.def("bins", &get_bins_from_file<hictk::cooler::File>);
  cooler_file.def("attributes", &cooler::get_cooler_attrs);

  cooler_file.def("fetch", &cooler::fetch, py::arg("range1") = "", py::arg("range2") = "",
                  py::arg("normalization") = "NONE", py::arg("count_type") = "int",
                  py::arg("join") = false, py::arg("query_type") = "UCSC");
  cooler_file.def("fetch_sparse", &cooler::fetch_sparse, py::arg("range1") = "",
                  py::arg("range2") = "", py::arg("normalization") = "NONE",
                  py::arg("count_type") = "int", py::arg("query_type") = "UCSC");

  return cooler;
}

static pybind11::module_ declare_hic_submodule(pybind11::module_& m) {
  auto hic = m.def_submodule("hic");
  auto hic_utils = hic.def_submodule("utils");

  hic_utils.def("is_hic_file", &hictk::hic::utils::is_hic_file,
                "test whether path points to a .hic file");

  auto hic_file = py::class_<hictk::hic::File>(hic, "File")
                      .def(py::init(&hic::file_ctor), py::arg("path"), py::arg("resolution"),
                           py::arg("matrix_type") = "observed", py::arg("matrix_unit") = "BP");

  hic_file.def("path", &hictk::hic::File::url);
  hic_file.def("name", &hictk::hic::File::name);
  hic_file.def("version", &hictk::hic::File::version);

  hic_file.def("bin_size", &hictk::hic::File::resolution);
  hic_file.def("nbins", &hictk::hic::File::nbins);
  hic_file.def("nchroms", &hictk::hic::File::nchroms);

  hic_file.def("chromosomes", &get_chromosomes_from_file<hictk::hic::File>);
  hic_file.def("bins", &get_bins_from_file<hictk::hic::File>);

  hic_file.def("fetch", &hic::fetch, py::arg("range1") = "", py::arg("range2") = "",
               py::arg("normalization") = "NONE", py::arg("count_type") = "int",
               py::arg("join") = false, py::arg("query_type") = "UCSC");
  hic_file.def("fetch_sparse", &hic::fetch_sparse, py::arg("range1") = "", py::arg("range2") = "",
               py::arg("normalization") = "NONE", py::arg("count_type") = "int",
               py::arg("query_type") = "UCSC");

  return hic;
}

static void declare_file_class(pybind11::module_& m) {
  auto file = py::class_<hictk::File>(m, "File").def(
      py::init(&file::ctor), py::arg("path"), py::arg("resolution"),
      py::arg("matrix_type") = "observed", py::arg("matrix_unit") = "BP");

  file.def("uri", &hictk::File::uri);
  file.def("path", &hictk::File::path);

  file.def("is_hic", &hictk::File::is_hic);
  file.def("is_cooler", &hictk::File::is_cooler);

  file.def("chromosomes", &get_chromosomes_from_file<hictk::File>);
  file.def("bins", &get_bins_from_file<hictk::File>);

  file.def("bin_size", &hictk::File::bin_size);
  file.def("nbins", &hictk::File::nbins);
  file.def("nchroms", &hictk::File::nchroms);

  file.def("fetch", &file::fetch, py::arg("range1") = "", py::arg("range2") = "",
           py::arg("normalization") = "NONE", py::arg("count_type") = "int",
           py::arg("join") = false, py::arg("query_type") = "UCSC");
  file.def("fetch_sparse", &file::fetch_sparse, py::arg("range1") = "", py::arg("range2") = "",
           py::arg("normalization") = "NONE", py::arg("count_type") = "int",
           py::arg("query_type") = "UCSC");
}

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(hictkpy, m) {
  [[maybe_unused]] auto np = py::module::import("numpy");
  [[maybe_unused]] auto pd = py::module::import("pandas");
  [[maybe_unused]] auto ss = py::module::import("scipy.sparse");
  m.attr("__version__") = hictk::config::version::str();

  m.doc() = "Blazing fast toolkit to work with .hic and .cool files";

  declare_cooler_submodule(m);
  declare_hic_submodule(m);
  declare_file_class(m);
}

}  // namespace hictkpy
