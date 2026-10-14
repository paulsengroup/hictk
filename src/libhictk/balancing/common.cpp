// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/common.hpp"

#include <fmt/format.h>

#include <stdexcept>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"

namespace hictk::balancing::internal {

void check_storage_mode(const cooler::File& clr) {
  if (clr.attributes().storage_mode != "symmetric-upper") {
    throw std::runtime_error(fmt::format(
        FMT_STRING("balancing interactions from files with storage-mode=\"{}\" is not supported"),
        clr.attributes().storage_mode.value_or("unknown")));
  }
}

void check_storage_mode(const File& f) {
  std::visit([](const auto& ff) { check_storage_mode(ff); }, f.get());
}

void check_bin_type(const BinTable& bins) {
  if (bins.type() == BinTable::Type::variable) {
    throw std::runtime_error(
        "balancing interactions from files with variable bin sizes is not supported");
  }
}

}  // namespace hictk::balancing::internal
