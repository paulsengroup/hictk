// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
// clang-format on

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_NULL_DEREF
#include <highfive/H5File.hpp>
HICTK_DISABLE_WARNING_POP
#include <highfive/H5Group.hpp>
#include <string>

namespace hictk::cooler {

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
struct RootGroup {
  HighFive::Group group{};

  [[nodiscard]] constexpr HighFive::Group &operator()() noexcept { return group; };
  [[nodiscard]] constexpr const HighFive::Group &operator()() const noexcept { return group; };

  [[nodiscard]] inline std::string file_name() const { return group.getFile().getName(); }
  [[nodiscard]] inline std::string hdf5_path() const { return group.getPath(); }
  [[nodiscard]] inline std::string uri() const {
    return fmt::format(FMT_STRING("{}::{}"), file_name(), hdf5_path());
  }
};

struct Group {
  RootGroup root_group{};
  HighFive::Group group{};

  [[nodiscard]] constexpr HighFive::Group &operator()() noexcept { return group; };
  [[nodiscard]] constexpr const HighFive::Group &operator()() const noexcept { return group; };
};
HICTK_DISABLE_WARNING_POP

using GroupMap = phmap::flat_hash_map<std::string, Group>;

}  // namespace hictk::cooler
