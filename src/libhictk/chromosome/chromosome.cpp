// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/chromosome.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/hash.hpp"

namespace hictk {

Chromosome::Chromosome(std::uint32_t id, std::string name_, std::uint32_t size_) noexcept
    : Chromosome(id, std::make_shared<const std::string>(std::move(name_)), size_) {}

Chromosome::Chromosome(std::uint32_t id, std::shared_ptr<const std::string> name_,
                       std::uint32_t size_)
    : _name(!!name_ ? std::move(name_) : std::make_shared<const std::string>("")),
      _id(id),
      _size(size_) {
  assert(_id != (std::numeric_limits<std::uint32_t>::max)());
  assert(_size != 0);
}

std::string_view Chromosome::name() const noexcept {
  return !!name_ptr() ? std::string_view{*name_ptr()} : "";  // NOLINT
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
bool Chromosome::is_all() const noexcept {
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

bool Chromosome::operator==(const Chromosome& other) const noexcept {
  return id() == other.id() && name() == other.name() && size() == other.size();
}

bool Chromosome::operator!=(const Chromosome& other) const noexcept { return !(*this == other); }

bool operator==(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() == b_name;
}
bool operator!=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() != b_name;
}

bool operator==(std::string_view a_name, const Chromosome& b) noexcept { return b == a_name; }
bool operator!=(std::string_view a_name, const Chromosome& b) noexcept { return !(b == a_name); }

}  // namespace hictk

std::size_t std::hash<hictk::Chromosome>::operator()(const hictk::Chromosome& c) const noexcept {
  return hictk::internal::hash_combine(0, c.id(), c.name(), c.size());
}
