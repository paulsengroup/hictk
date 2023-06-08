// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <functional>
#include <iterator>

namespace hicxx::internal {

constexpr HiCHeader::operator bool() const noexcept { return masterIndexOffset >= 0; }

inline bool HiCHeader::operator==(const HiCHeader &other) const noexcept {
    return url == other.url && masterIndexOffset == other.masterIndexOffset;
}

inline bool HiCHeader::operator!=(const HiCHeader &other) const noexcept {
    return !(*this == other);
}

inline std::size_t HiCHeader::nChromosomes() const noexcept { return chromosomes.size(); }

inline std::size_t HiCHeader::nResolutions() const noexcept { return resolutions.size(); }

inline const chromosome &HiCHeader::getChromosome(std::int32_t id) const noexcept {
    // Chromosomes are sorted by id, so we can use simple arithmetic on iterators to find the
    // chromosome with the given id
    assert(id < static_cast<std::int32_t>(chromosomes.size()));

    const auto it = std::next(chromosomes.begin(), id);
    assert(it->second.index == id);
    return it->second;
}

}  // namespace hicxx::internal

template <>
struct std::hash<hicxx::internal::HiCHeader> {
    inline std::size_t operator()(hicxx::internal::HiCHeader const &h) const noexcept {
        return hicxx::internal::hash_combine(0, h.url, h.masterIndexOffset);
    }
};
