// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/fuzzer/cooler.hpp"

#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>

namespace hictk::fuzzer::cooler {

namespace py = pybind11;

static py::module_ import_cooler() {
  try {
    return py::module_::import("cooler");
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("unable to import cooler: {}"), e.what()));
  }
}

[[nodiscard]] static py::object open_cooler(std::string_view uri) {
  try {
    return import_cooler().attr("Cooler")(uri);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to open Cooler at {}: {}"), uri, e.what()));
  }
}

std::string_view version() {
  auto importlib_metadata = py::module_::import("importlib.metadata");
  return py::cast<std::string_view>(importlib_metadata.attr("version")("cooler"));
}

Cooler::Cooler(std::string_view uri) : _clr(open_cooler(uri)) {}

std::string Cooler::uri() const noexcept {
  if (_clr) {
    return _clr.attr("uri").cast<std::string>();
  }
  return "";
}

}  // namespace hictk::fuzzer::cooler
