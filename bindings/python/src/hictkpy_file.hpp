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

inline py::object fetch(const hictk::File& f, std::string_view range1, std::string_view range2,
                        std::string_view normalization, std::string_view count_type, bool join,
                        std::string_view query_type) {
  return std::visit(
      [&](const auto& ff) -> py::object {
        return file_fetch(ff, range1, range2, normalization, count_type, join, query_type);
      },
      f.get());
}

}  // namespace hictkpy::file
