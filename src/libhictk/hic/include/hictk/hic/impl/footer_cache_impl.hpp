// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

#include "hictk/hic/footer.hpp"

namespace hictk::hic::internal {

inline bool FooterCache::HiCFooterPtrCmp::operator()(
    const std::shared_ptr<const HiCFooter> &f1,
    const std::shared_ptr<const HiCFooter> &f2) const noexcept {
  if (!f1 || !f2) {
    return f1 == f2;
  }
  return *f1 == *f2;
}

inline bool FooterCache::HiCFooterPtrCmp::operator()(
    const HiCFooterMetadata &m1, const std::shared_ptr<const HiCFooter> &f2) const noexcept {
  if (!f2) {
    return false;
  }
  return m1 == f2->metadata();
}

inline bool FooterCache::HiCFooterPtrCmp::operator()(const std::shared_ptr<const HiCFooter> &f1,
                                                     const HiCFooterMetadata &m2) const noexcept {
  if (!f1) {
    return false;
  }
  return f1->metadata() == m2;
}

inline std::size_t FooterCache::HiCFooterPtrHasher::operator()(
    const std::shared_ptr<const HiCFooter> &f) const noexcept {
  if (!f) {
    return std::hash<const HiCFooter *>{}(f.get());
  }
  return (*this)(f->metadata());
}
inline std::size_t FooterCache::HiCFooterPtrHasher::operator()(
    const HiCFooterMetadata &m) const noexcept {
  return std::hash<HiCFooterMetadata>{}(m);
}

inline auto FooterCache::begin() const noexcept -> decltype(_footers.cbegin()) {
  return _footers.begin();
}
inline auto FooterCache::end() const noexcept -> decltype(_footers.cbegin()) {
  return _footers.end();
}

inline auto FooterCache::cbegin() const noexcept -> decltype(_footers.cbegin()) {
  return _footers.cbegin();
}
inline auto FooterCache::cend() const noexcept -> decltype(_footers.cbegin()) {
  return _footers.cend();
}

inline auto FooterCache::emplace(HiCFooter f) -> decltype(_footers.emplace()) {
  return _footers.emplace(std::make_shared<const HiCFooter>(std::move(f)));
}
inline auto FooterCache::find(const HiCFooterMetadata &m) -> const_iterator {
  return _footers.find(m);
}
inline std::size_t FooterCache::size() const noexcept { return _footers.size(); }
inline void FooterCache::clear() { _footers.clear(); }

}  // namespace hictk::hic::internal
