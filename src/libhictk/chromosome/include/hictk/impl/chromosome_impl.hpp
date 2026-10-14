// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace hictk {

constexpr Chromosome::operator bool() const noexcept { return id() != null_id; }

constexpr std::uint32_t Chromosome::id() const noexcept { return _id; }

constexpr const std::shared_ptr<const std::string>& Chromosome::name_ptr() const noexcept {
  return _name;
}

constexpr std::uint32_t Chromosome::size() const noexcept { return _size; }

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
