// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

#include <string>
#include <string_view>
#include <utility>

namespace hictk::cooler {

struct CoolerURI {
  std::string file_path;
  std::string group_path;

  CoolerURI() = default;
  CoolerURI(std::string_view p1, std::string_view p2);
  CoolerURI(std::string p1, std::string p2);
  explicit CoolerURI(std::pair<std::string_view, std::string_view> paths);
  explicit CoolerURI(std::pair<std::string, std::string> paths);
  // clang-format on
};

[[nodiscard]] CoolerURI parse_cooler_uri(std::string_view uri);
}  // namespace hictk::cooler

#include "./impl/uri_impl.hpp"  // NOLINT
