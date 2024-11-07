// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <cstdint>
#include <string>
#include <vector>

#include "hictk/binary_buffer.hpp"
#include "hictk/filestream.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

struct HiCHeader {
  std::string url{};
  std::int32_t version{-1};
  std::int64_t footerPosition{-1};
  std::string genomeID{};
  std::int64_t normVectorIndexPosition{-1};
  std::int64_t normVectorIndexLength{-1};
  Reference chromosomes{};
  std::vector<std::uint32_t> resolutions{};
  phmap::flat_hash_map<std::string, std::string> attributes{};

  constexpr explicit operator bool() const noexcept;
  bool operator==(const HiCHeader& other) const noexcept;
  bool operator!=(const HiCHeader& other) const noexcept;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static HiCHeader deserialize(std::streampos offset, filestream::FileStream<>& fs);
  [[nodiscard]] static HiCHeader unsafe_deserialize(std::streampos offset,
                                                    filestream::FileStream<>& fs);
};

}  // namespace hictk::hic::internal

#include "./impl/header_impl.hpp"  // NOLINT
