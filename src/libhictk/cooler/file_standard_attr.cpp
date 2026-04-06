// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include <fmt/format.h>
#include <fmt/chrono.h>

#include "hictk/cooler/cooler.hpp"
// clang-format on

#include <cassert>
#include <cstdint>
#include <ctime>
#include <string>

#include "hictk/bin_table.hpp"

namespace hictk::cooler {

Attributes Attributes::init_empty() noexcept {
  Attributes attrs{};

  attrs.bin_type = BinTable::Type::fixed;
  attrs.creation_date.reset();
  attrs.format_url.reset();
  attrs.generated_by.reset();
  attrs.assembly.reset();
  attrs.nbins.reset();
  attrs.nchroms.reset();
  attrs.metadata.reset();
  attrs.storage_mode.reset();
  attrs.sum.reset();
  attrs.cis.reset();

  return attrs;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
bool Attributes::operator==(const Attributes& other) const noexcept {
  if (!sum.has_value()) {
    assert(!sum->valueless_by_exception());  // NOLINT(*-unchecked-optional-access)
  }
  if (!other.sum.has_value()) {
    assert(!other.sum->valueless_by_exception());  // NOLINT(*-unchecked-optional-access)
  }
  // clang-format off
  return bin_size == other.bin_size &&
         bin_type == other.bin_type &&
         format == other.format &&
         format_version == other.format_version &&
         storage_mode == other.storage_mode &&
         creation_date == other.creation_date &&
         generated_by == other.generated_by &&
         assembly == other.assembly &&
         metadata == other.metadata &&
         format_url == other.format_url &&
         nbins == other.nbins &&
         nchroms == other.nchroms &&
         nnz == other.nnz &&
         sum == other.sum &&
         cis == other.cis;
  // clang-format on
}

bool Attributes::operator!=(const Attributes& other) const noexcept { return !(*this == other); }

std::string Attributes::generate_creation_date() {
  // e.g. 2022-07-26T20:35:19
  return fmt::format(FMT_STRING("{:%FT%T}"), fmt::gmtime(std::time(nullptr)));
}
}  // namespace hictk::cooler
