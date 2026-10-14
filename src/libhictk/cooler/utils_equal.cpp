// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>

#include "hictk/cooler/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/utils.hpp"

namespace {
bool attributes_are_equal(hictk::cooler::Attributes attr1, hictk::cooler::Attributes attr2) {
  attr1.creation_date = "";  // NOLINT
  attr2.creation_date = "";  // NOLINT
  attr1.metadata = "";       // NOLINT
  attr2.metadata = "";       // NOLINT

  return attr1 == attr2;
}

template <typename T = std::int64_t>
bool datasets_are_equal(const hictk::cooler::Dataset& d1, const hictk::cooler::Dataset& d2) {
  if (d1.size() != d2.size()) {
    return false;
  }
  if (d1.empty()) {
    return true;
  }

  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  return std::equal(d1.begin<T>(256'000), d1.end<T>(256'000), d2.begin<T>(256'000));
}

}  // namespace

namespace hictk::cooler::utils {

bool equal(std::string_view uri1, std::string_view uri2, bool ignore_attributes) {
  if (uri1 == uri2) {
    return true;
  }
  return equal(File::open_read_once(uri1), File::open_read_once(uri2), ignore_attributes);
}

bool equal(const File& clr1, const File& clr2, bool ignore_attributes) {
  if (clr1.uri() == clr2.uri()) {
    return true;
  }

  if (!ignore_attributes && !attributes_are_equal(clr1.attributes(), clr2.attributes())) {
    return false;
  }

  const auto float_counts = clr1.has_float_pixels() || clr2.has_float_pixels();
  for (const auto name : MANDATORY_DATASET_NAMES) {
    bool difference_found;  // NOLINT
    if (name == "chroms/name") {
      difference_found = !datasets_are_equal<std::string>(clr1.dataset(name), clr2.dataset(name));
    } else if (name == "pixels/count" && float_counts) {
      difference_found = !datasets_are_equal<double>(clr1.dataset(name), clr2.dataset(name));

    } else {
      difference_found = !datasets_are_equal(clr1.dataset(name), clr2.dataset(name));
    }

    if (difference_found) {
      // fmt::print(FMT_STRING("Difference found in {} ({}, {})!\n"), name, clr1.uri(), clr2.uri());
      return false;
    }
  }
  return true;
}

}  // namespace hictk::cooler::utils
