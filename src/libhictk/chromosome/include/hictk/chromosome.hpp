// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace hictk {

class Chromosome {
  static constexpr std::uint32_t null_id{(std::numeric_limits<std::uint32_t>::max)()};

  std::shared_ptr<std::string> _name{};
  std::uint32_t _id{null_id};
  std::uint32_t _size{};

 public:
  Chromosome() = default;
  Chromosome(std::uint32_t id_, std::string name_, std::uint32_t size_) noexcept;

  [[nodiscard]] constexpr explicit operator bool() const noexcept;

  [[nodiscard]] constexpr std::uint32_t id() const noexcept;
  [[nodiscard]] std::string_view name() const noexcept;
  [[nodiscard]] constexpr std::uint32_t size() const noexcept;
  [[nodiscard]] bool is_all() const noexcept;

  [[nodiscard]] constexpr bool operator<(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator==(const Chromosome& other) const noexcept;
  [[nodiscard]] bool operator!=(const Chromosome& other) const noexcept;

  friend bool operator==(const Chromosome& a, std::string_view b_name) noexcept;
  friend bool operator!=(const Chromosome& a, std::string_view b_name) noexcept;

  friend bool operator==(std::string_view a_name, const Chromosome& b) noexcept;
  friend bool operator!=(std::string_view a_name, const Chromosome& b) noexcept;

  friend constexpr bool operator<(const Chromosome& a, std::uint32_t b_id) noexcept;
  friend constexpr bool operator>(const Chromosome& a, std::uint32_t b_id) noexcept;
  friend constexpr bool operator<=(const Chromosome& a, std::uint32_t b_id) noexcept;
  friend constexpr bool operator>=(const Chromosome& a, std::uint32_t b_id) noexcept;
  friend constexpr bool operator==(const Chromosome& a, std::uint32_t b_id) noexcept;
  friend constexpr bool operator!=(const Chromosome& a, std::uint32_t b_id) noexcept;

  friend constexpr bool operator<(std::uint32_t a_id, const Chromosome& b) noexcept;
  friend constexpr bool operator>(std::uint32_t a_id, const Chromosome& b) noexcept;
  friend constexpr bool operator<=(std::uint32_t a_id, const Chromosome& b) noexcept;
  friend constexpr bool operator>=(std::uint32_t a_id, const Chromosome& b) noexcept;
  friend constexpr bool operator==(std::uint32_t a_id, const Chromosome& b) noexcept;
  friend constexpr bool operator!=(std::uint32_t a_id, const Chromosome& b) noexcept;
};

struct ChromosomeCmp {
  using is_transparent = void;

  constexpr bool operator()(const Chromosome& c1, const Chromosome& c2) const noexcept;
  constexpr bool operator()(std::uint32_t id1, const Chromosome& c2) const noexcept;
  constexpr bool operator()(const Chromosome& c1, std::uint32_t id2) const noexcept;
};

}  // namespace hictk

namespace std {
template <>
struct hash<hictk::Chromosome> {
  size_t operator()(const hictk::Chromosome& c) const;
};
}  // namespace std

#include "../../chromosome_impl.hpp"
