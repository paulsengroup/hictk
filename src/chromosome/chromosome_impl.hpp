// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/common.hpp"
#include "hictk/hash.hpp"

namespace hictk {

inline Chromosome::Chromosome(std::uint32_t id, std::string name_, std::uint32_t size_) noexcept
    : _name(std::move(name_)), _id(id), _size(size_) {
  assert(_id != (std::numeric_limits<std::uint32_t>::max)());
  assert(_size != 0);
}

constexpr Chromosome::operator bool() const noexcept { return this->id() != Chromosome::null_id; }

constexpr std::uint32_t Chromosome::id() const noexcept { return this->_id; }

inline std::string_view Chromosome::name() const noexcept { return this->_name; }

constexpr std::uint32_t Chromosome::size() const noexcept { return this->_size; }

inline bool Chromosome::is_all() const noexcept {
  if (this->name() == "All") {
    return true;
  }

  std::string name{this->name()};

  std::transform(name.begin(), name.end(), name.begin(),
                 [&](const char c) { return std::tolower(c); });

  return name == "all";
}

constexpr bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return this->id() < other.id();
}

constexpr bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return this->id() > other.id();
}

constexpr bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return this->id() <= other.id();
}

constexpr bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return this->id() >= other.id();
}

inline bool Chromosome::operator==(const Chromosome& other) const noexcept {
  return this->id() == other.id() && this->name() == other.name() && this->size() == other.size();
}

inline bool Chromosome::operator!=(const Chromosome& other) const noexcept {
  return !(*this == other);
}

inline bool operator==(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() == b_name;
}
inline bool operator!=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() != b_name;
}

inline bool operator==(std::string_view a_name, const Chromosome& b) noexcept {
  return b == a_name;
}
inline bool operator!=(std::string_view a_name, const Chromosome& b) noexcept {
  return !(b == a_name);
}

constexpr bool operator<(const Chromosome& a, std::uint32_t b_id) noexcept { return a.id() < b_id; }
constexpr bool operator>(const Chromosome& a, std::uint32_t b_id) noexcept { return a.id() > b_id; }
constexpr bool operator<=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() <= b_id;
}
constexpr bool operator>=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() >= b_id;
}
constexpr bool operator==(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() == b_id;
}
constexpr bool operator!=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() != b_id;
}

constexpr bool operator<(std::uint32_t a_id, const Chromosome& b) noexcept { return b > a_id; }
constexpr bool operator>(std::uint32_t a_id, const Chromosome& b) noexcept { return b < a_id; }
constexpr bool operator<=(std::uint32_t a_id, const Chromosome& b) noexcept { return b >= a_id; }
constexpr bool operator>=(std::uint32_t a_id, const Chromosome& b) noexcept { return b <= a_id; }
constexpr bool operator==(std::uint32_t a_id, const Chromosome& b) noexcept { return b == a_id; }
constexpr bool operator!=(std::uint32_t a_id, const Chromosome& b) noexcept { return b != a_id; }

}  // namespace hictk

inline std::size_t std::hash<hictk::Chromosome>::operator()(const hictk::Chromosome& c) const {
  return hictk::internal::hash_combine(0, c.id(), c.name(), c.size());
}
