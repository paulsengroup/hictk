// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/chromosome.hpp"
#include "hictk/hash.hpp"

namespace hictk {

inline Chromosome::Chromosome(std::uint32_t id, std::string name_, std::uint32_t size_) noexcept
    : _name(std::make_shared<std::string>(std::move(name_))), _id(id), _size(size_) {
  assert(_id != (std::numeric_limits<std::uint32_t>::max)());
  assert(_size != 0);
}

constexpr Chromosome::operator bool() const noexcept { return id() != Chromosome::null_id; }

constexpr std::uint32_t Chromosome::id() const noexcept { return _id; }

inline std::string_view Chromosome::name() const noexcept {
  return !!_name ? std::string_view{*_name} : "";  // NOLINT
}

constexpr std::uint32_t Chromosome::size() const noexcept { return _size; }

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
inline bool Chromosome::is_all() const noexcept {
  constexpr std::string_view all{"All"};
  if (name() == all) {
    return true;
  }

  if (name().size() != all.size()) {
    return false;
  }

  return std::equal(name().begin(), name().end(), all.begin(),
                    [](const char a, const char b) { return std::tolower(a) == std::tolower(b); });
}

constexpr bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return id() < other.id();
}

constexpr bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return id() > other.id();
}

constexpr bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return id() <= other.id();
}

constexpr bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return id() >= other.id();
}

inline bool Chromosome::operator==(const Chromosome& other) const noexcept {
  return id() == other.id() && name() == other.name() && size() == other.size();
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

constexpr bool ChromosomeCmp::operator()(const Chromosome& c1,
                                         const Chromosome& c2) const noexcept {
  return c1 < c2;
}
constexpr bool ChromosomeCmp::operator()(std::uint32_t id1, const Chromosome& c2) const noexcept {
  return id1 < c2.id();
}
constexpr bool ChromosomeCmp::operator()(const Chromosome& c1, std::uint32_t id2) const noexcept {
  return c1.id() < id2;
}

}  // namespace hictk

inline std::size_t std::hash<hictk::Chromosome>::operator()(const hictk::Chromosome& c) const {
  return hictk::internal::hash_combine(0, c.id(), c.name(), c.size());
}
