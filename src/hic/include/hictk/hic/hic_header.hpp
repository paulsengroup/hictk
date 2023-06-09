// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "hictk/hic/common.hpp"
#include "hictk/reference.hpp"

namespace hictk::internal {

struct HiCHeader {
  std::string url{};
  std::int32_t version{-1};
  std::int64_t masterIndexOffset{-1};
  std::string genomeID{};
  std::int64_t nviPosition{-1};
  std::int64_t nviLength{-1};
  Reference chromosomes{};
  std::vector<std::uint32_t> resolutions{};

  constexpr explicit operator bool() const noexcept;
  bool operator==(const HiCHeader &other) const noexcept;
  bool operator!=(const HiCHeader &other) const noexcept;
};

}  // namespace hictk::internal

#include "../../../hic_header_impl.hpp"
