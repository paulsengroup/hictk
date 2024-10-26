// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_NULL_DEREFERENCE
#include <highfive/H5File.hpp>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>

#include <highfive/H5Group.hpp>
#include <optional>
#include <string>

namespace hictk::cooler {

class RootGroup {
  std::optional<HighFive::Group> _group{};

 public:
  RootGroup() = default;
  explicit RootGroup(HighFive::Group grp) noexcept : _group(std::move(grp)) {}

  RootGroup(const RootGroup& other) = default;
  // NOLINTNEXTLINE(bugprone-exception-escape)
  RootGroup(RootGroup&& other) noexcept : _group(std::move(other._group)) {}

  ~RootGroup() noexcept = default;

  RootGroup& operator=(const RootGroup& other) = default;
  // NOLINTNEXTLINE(bugprone-exception-escape)
  RootGroup& operator=(RootGroup&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    _group = std::move(other._group);

    return *this;
  }

  [[nodiscard]] constexpr HighFive::Group& operator()() {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return _group.value();
  }
  [[nodiscard]] constexpr const HighFive::Group& operator()() const {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return _group.value();
  }

  [[nodiscard]] std::string file_name() const {
    if (!_group.has_value()) {
      return "";
    }
    return _group->getFile().getName();
  }
  [[nodiscard]] std::string hdf5_path() const {
    if (!_group.has_value()) {
      return "";
    }
    return _group->getPath();
  }
  [[nodiscard]] std::string uri() const {
    if (!_group.has_value()) {
      return "";
    }
    return fmt::format(FMT_STRING("{}::{}"), file_name(), hdf5_path());
  }
};

class Group {
  RootGroup _root_group{};
  std::optional<HighFive::Group> _group{};

 public:
  Group() = default;
  Group(RootGroup root_grp, HighFive::Group grp) noexcept
      : _root_group(std::move(root_grp)), _group(std::move(grp)) {}
  [[nodiscard]] constexpr HighFive::Group& operator()() {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return _group.value();
  }
  [[nodiscard]] constexpr const HighFive::Group& operator()() const {
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return _group.value();
  }
};

using GroupMap = phmap::flat_hash_map<std::string, Group>;

}  // namespace hictk::cooler
