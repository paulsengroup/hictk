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

#include <cstddef>
#include <memory>

#include "hictk/hic/footer.hpp"

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

  using FooterMap =
      phmap::flat_hash_set<std::shared_ptr<const HiCFooter>, HiCFooterPtrHasher, HiCFooterPtrCmp>;
  FooterMap _footers;

 public:
  using difference_type = FooterMap::difference_type;
  using iterator = FooterMap::iterator;
  using const_iterator = FooterMap::iterator;
  FooterCache() = default;

  [[nodiscard]] auto begin() const noexcept -> decltype(_footers.cbegin());
  [[nodiscard]] auto end() const noexcept -> decltype(_footers.cbegin());

  [[nodiscard]] auto cbegin() const noexcept -> decltype(_footers.cbegin());
  [[nodiscard]] auto cend() const noexcept -> decltype(_footers.cbegin());

  [[nodiscard]] auto emplace(HiCFooter f) -> decltype(_footers.emplace());
  [[nodiscard]] auto find(const HiCFooterMetadata& m) -> const_iterator;

  [[nodiscard]] std::size_t size() const noexcept;
  void clear();
};

}  // namespace hictk::hic::internal

#include "./impl/footer_cache_impl.hpp"  // NOLINT
