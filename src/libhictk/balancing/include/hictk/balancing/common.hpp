// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <stdexcept>
#include <type_traits>
#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"

namespace hictk::balancing::internal {
inline void check_storage_mode(const cooler::File& clr) {
  if (clr.attributes().storage_mode != "symmetric-upper") {
    throw std::runtime_error(fmt::format(
        FMT_STRING("balancing interactions from files with storage-mode=\"{}\" is not supported"),
        clr.attributes().storage_mode.value_or("unknown")));
  }
}

template <typename File>
inline void check_storage_mode([[maybe_unused]] const File& f) {
  if constexpr (std::is_same_v<File, cooler::File>) {
    check_storage_mode(f);
    return;
  }

  if constexpr (std::is_same_v<File, hictk::File>) {
    if (std::holds_alternative<cooler::File>(f.get())) {
      check_storage_mode(std::get<cooler::File>(f.get()));
    }
  }
}

}  // namespace hictk::balancing::internal
