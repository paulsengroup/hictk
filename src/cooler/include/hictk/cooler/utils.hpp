// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/cooler.hpp"

namespace hictk::cooler::utils {

/// Iterable of hictk::File or strings
template <typename Str>
void merge(Str first_file, Str last_file, std::string_view dest_uri,
           bool overwrite_if_exists = false, std::size_t chunk_size = 500'000, bool quiet = true);

[[nodiscard]] bool equal(std::string_view uri1, std::string_view uri2,
                         bool ignore_attributes = true);
[[nodiscard]] bool equal(const File& clr1, const File& clr2, bool ignore_attributes = true);

}  // namespace hictk::cooler::utils

#include "../../../utils_equal_impl.hpp"
#include "../../../utils_merge_impl.hpp"
