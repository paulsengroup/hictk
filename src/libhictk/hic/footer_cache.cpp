// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/footer_cache.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/hic/footer.hpp"

namespace hictk::hic::internal {

bool FooterCache::HiCFooterPtrCmp::operator()(
    const std::shared_ptr<const HiCFooter> &f1,
    const std::shared_ptr<const HiCFooter> &f2) const noexcept {
  if (!f1 || !f2) {
    return f1 == f2;
  }
  return *f1 == *f2;
}

bool FooterCache::HiCFooterPtrCmp::operator()(
    const HiCFooterMetadata &m1, const std::shared_ptr<const HiCFooter> &f2) const noexcept {
  if (!f2) {
    return false;
  }
  return m1 == f2->metadata();
}

bool FooterCache::HiCFooterPtrCmp::operator()(const std::shared_ptr<const HiCFooter> &f1,
                                              const HiCFooterMetadata &m2) const noexcept {
  if (!f1) {
    return false;
  }
  return f1->metadata() == m2;
}

std::size_t FooterCache::HiCFooterPtrHasher::operator()(
    const std::shared_ptr<const HiCFooter> &f) const noexcept {
  if (!f) {
    return std::hash<const HiCFooter *>{}(f.get());
  }
  return (*this)(f->metadata());
}
std::size_t FooterCache::HiCFooterPtrHasher::operator()(const HiCFooterMetadata &m) const noexcept {
  return std::hash<HiCFooterMetadata>{}(m);
}

auto FooterCache::begin() const noexcept -> decltype(_footers.cbegin()) { return _footers.begin(); }
auto FooterCache::end() const noexcept -> decltype(_footers.cbegin()) { return _footers.end(); }

auto FooterCache::cbegin() const noexcept -> decltype(_footers.cbegin()) {
  return _footers.cbegin();
}
auto FooterCache::cend() const noexcept -> decltype(_footers.cbegin()) { return _footers.cend(); }

auto FooterCache::emplace(HiCFooter f) -> decltype(_footers.emplace()) {
  return _footers.emplace(std::make_shared<const HiCFooter>(std::move(f)));
}
auto FooterCache::find(const HiCFooterMetadata &m) -> const_iterator { return _footers.find(m); }
std::size_t FooterCache::size() const noexcept { return _footers.size(); }
void FooterCache::clear() { _footers.clear(); }

}  // namespace hictk::hic::internal
