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
  explicit RootGroup(HighFive::Group grp) noexcept;

  RootGroup(const RootGroup& other) = default;
  // NOLINTNEXTLINE(bugprone-exception-escape)
  RootGroup(RootGroup&& other) noexcept;

  ~RootGroup() noexcept = default;

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] constexpr bool operator!() const noexcept;

  RootGroup& operator=(const RootGroup& other) = default;
  RootGroup& operator=(RootGroup&& other) noexcept;

  [[nodiscard]] constexpr HighFive::Group& operator()();
  [[nodiscard]] constexpr const HighFive::Group& operator()() const;

  [[nodiscard]] std::string file_name() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string uri() const;

  template <typename N>
  void read_attribute(std::string_view key, N& buff) const;
  template <typename N>
  [[nodiscard]] N read_attribute(std::string_view key) const;

  template <typename N>
  [[nodiscard]] bool try_read_attribute(std::string_view key, N& buff) const;
  template <typename N>
  [[nodiscard]] std::optional<N> try_read_attribute(std::string_view key) const;

  template <typename N>
  [[nodiscard]] bool try_read_sum_attribute(std::string_view key, N& buff) const;
  template <typename N>
  [[nodiscard]] std::optional<N> try_read_sum_attribute(std::string_view key) const;
};

class Group {
  RootGroup _root_group{};
  std::optional<HighFive::Group> _group{};

 public:
  Group() = default;
  Group(RootGroup root_grp, HighFive::Group grp) noexcept;
  [[nodiscard]] constexpr HighFive::Group& operator()();
  [[nodiscard]] constexpr const HighFive::Group& operator()() const;
};

using GroupMap = phmap::flat_hash_map<std::string, Group>;

}  // namespace hictk::cooler

#include "./impl/group_impl.hpp"
