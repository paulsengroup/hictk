// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "./common.hpp"
#include "hictk/file.hpp"

namespace hictkpy::file {
[[nodiscard]] hictk::File ctor(std::string_view path, std::int32_t resolution,
                               std::string_view matrix_type, std::string_view matrix_unit) {
  return hictk::File{std::string{path}, static_cast<std::uint32_t>(resolution),
                     hictk::hic::ParseMatrixTypeStr(std::string{matrix_type}),
                     hictk::hic::ParseUnitStr(std::string{matrix_unit})};
}

static py::object fetch(const hictk::File& f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, bool join, std::string_view query_type) {
  const auto qt =
      query_type == "UCSC" ? hictk::File::QUERY_TYPE::UCSC : hictk::File::QUERY_TYPE::BED;

  auto sel = range2.empty() || range1 == range2
                 ? f.fetch(range1, hictk::balancing::Method(normalization), qt)
                 : f.fetch(range1, range2, hictk::balancing::Method(normalization), qt);
  if (normalization == "NONE") {
    return pixel_iterators_to_df(f.bins(), sel.begin<std::int32_t>(), sel.end<std::int32_t>(),
                                 join);
  }
  return pixel_iterators_to_df(f.bins(), sel.begin<double>(), sel.end<double>(), join);
}
}  // namespace hictkpy::file
