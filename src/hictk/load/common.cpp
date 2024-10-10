// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./common.hpp"

#include <cassert>
#include <string_view>
#include <variant>

#include "hictk/type_traits.hpp"

namespace hictk::tools {

Stats& Stats::operator+=(const Stats& other) {
  std::visit(
      [&](auto& sum_) {
        using T = remove_cvref_t<decltype(sum_)>;

        sum_ += std::get<T>(other.sum);
      },
      sum);
  nnz += other.nnz;

  return *this;
}

Format format_from_string(std::string_view s) {
  if (s == "coo") {
    return Format::COO;
  }
  if (s == "bg2") {
    return Format::BG2;
  }
  if (s == "validpairs") {
    return Format::VP;
  }
  assert(s == "4dn");
  return Format::_4DN;
}

}  // namespace hictk::tools
