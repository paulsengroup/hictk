// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <parallel_hashmap/phmap.h>

#include <cstdint>
#include <string>
#include <vector>

#include "hictk/reference.hpp"

namespace hictk::hic::internal {

struct HiCHeader {
  std::string url{};
  std::int32_t version{-1};
  std::int64_t masterIndexOffset{-1};
  std::string genomeID{};
  std::int64_t nviPosition{-1};
  std::int64_t nviLength{-1};
  Reference chromosomes{};
  std::vector<std::uint32_t> resolutions{};
  phmap::flat_hash_map<std::string, std::string> attributes{};

  constexpr explicit operator bool() const noexcept;
  bool operator==(const HiCHeader &other) const noexcept;
  bool operator!=(const HiCHeader &other) const noexcept;
};

}  // namespace hictk::hic::internal

#include "./impl/header_impl.hpp"  // NOLINT
