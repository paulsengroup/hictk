// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <string_view>

#include "hictk/cooler.hpp"

namespace hictk::cooler::utils {

inline bool equal(std::string_view uri1, std::string_view uri2, bool ignore_attributes) {
  if (uri1 == uri2) {
    return true;
  }
  return equal(File::open_read_only_read_once(uri1), File::open_read_only_read_once(uri2),
               ignore_attributes);
}

namespace internal {
inline bool attributes_are_equal(StandardAttributes attr1, StandardAttributes attr2) {
  attr1.creation_date = "";  // NOLINT
  attr2.creation_date = "";  // NOLINT
  attr1.metadata = "";       // NOLINT
  attr2.metadata = "";       // NOLINT

  return attr1 == attr2;
}

template <typename T = std::int64_t>
inline bool datasets_are_equal(const Dataset& d1, const Dataset& d2) {
  if (d1.size() != d2.size()) {
    return false;
  }
  if (d1.empty()) {
    return true;
  }

  return std::equal(d1.begin<T>(256'000), d1.end<T>(256'000), d2.begin<T>(256'000));
}

}  // namespace internal

inline bool equal(const File& clr1, const File& clr2, bool ignore_attributes) {
  if (clr1.uri() == clr2.uri()) {
    return true;
  }

  if (!ignore_attributes && !internal::attributes_are_equal(clr1.attributes(), clr2.attributes())) {
    return false;
  }

  const auto float_counts = clr1.has_float_pixels() || clr2.has_float_pixels();
  for (const auto name : MANDATORY_DATASET_NAMES) {
    bool difference_found;  // NOLINT
    if (name == "chroms/name") {
      difference_found =
          !internal::datasets_are_equal<std::string>(clr1.dataset(name), clr2.dataset(name));
    } else if (name == "pixels/count" && float_counts) {
      difference_found =
          !internal::datasets_are_equal<double>(clr1.dataset(name), clr2.dataset(name));

    } else {
      difference_found = !internal::datasets_are_equal(clr1.dataset(name), clr2.dataset(name));
    }

    if (difference_found) {
      // fmt::print(FMT_STRING("Difference found in {} ({}, {})!\n"), name, clr1.uri(), clr2.uri());
      return false;
    }
  }
  return true;
}

}  // namespace hictk::cooler::utils
