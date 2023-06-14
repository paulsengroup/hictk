// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#include <tsl/ordered_map.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/hic_footer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class FooterCache {
  struct HiCFooterPtrCmp {
    using is_transparent = void;

    bool operator()(const std::shared_ptr<const HiCFooter>& f1,
                    const std::shared_ptr<const HiCFooter>& f2) const noexcept;
    bool operator()(const HiCFooterMetadata& m1,
                    const std::shared_ptr<const HiCFooter>& f2) const noexcept;
    bool operator()(const std::shared_ptr<const HiCFooter>& f1,
                    const HiCFooterMetadata& m2) const noexcept;
  };

  struct HiCFooterPtrHasher {
    using is_transparent = void;

    std::size_t operator()(const std::shared_ptr<const HiCFooter>& f) const noexcept;
    std::size_t operator()(const HiCFooterMetadata& m) const noexcept;
  };

  using MapT =
      phmap::flat_hash_set<std::shared_ptr<const HiCFooter>, HiCFooterPtrHasher, HiCFooterPtrCmp>;
  MapT _cache;

 public:
  using difference_type = MapT::difference_type;
  using iterator = MapT::iterator;
  using const_iterator = MapT::iterator;
  FooterCache() = default;

  auto begin() const noexcept -> decltype(_cache.cbegin());
  auto end() const noexcept -> decltype(_cache.cbegin());

  auto cbegin() const noexcept -> decltype(_cache.cbegin());
  auto cend() const noexcept -> decltype(_cache.cbegin());

  auto emplace(HiCFooter&& f) -> decltype(_cache.emplace());
  auto find(const HiCFooterMetadata& m) -> const_iterator;

  [[nodiscard]] std::size_t size() const noexcept;
  void clear();
};

}  // namespace hictk::hic::internal

#include "../../../footer_cache_impl.hpp"
