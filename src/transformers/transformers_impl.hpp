// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::transformers {

inline JoinGenomicCoords::JoinGenomicCoords(cooler::PixelSelector<> sel,
                                            std::shared_ptr<const BinTable> bins)
    : _sel(std::in_place_type<cooler::PixelSelector<>>, std::move(sel)), _bins(std::move(bins)) {}
inline JoinGenomicCoords::JoinGenomicCoords(hic::PixelSelector sel,
                                            std::shared_ptr<const BinTable> bins)
    : _sel(std::in_place_type<hic::PixelSelector>, std::move(sel)), _bins(std::move(bins)) {}

template <typename N>
inline auto JoinGenomicCoords::begin() -> iterator<N> {
  std::visit([&](auto& sel) { return iterator<N>{sel.template begin<N>(), _bins}; }, _sel);
}
template <typename N>
inline auto JoinGenomicCoords::end() -> iterator<N> {
  std::visit([&](auto& sel) { return iterator<N>{sel.template end<N>(), _bins}; }, _sel);
}

template <typename N>
inline auto JoinGenomicCoords::cbegin() -> iterator<N> {
  return this->begin<N>();
}
template <typename N>
inline auto JoinGenomicCoords::cend() -> iterator<N> {
  return this->end<N>();
}

template <typename N>
inline std::vector<Pixel<N>> JoinGenomicCoords::read_all() const {
  return std::visit([](const auto& sel) { return sel.template read_all<N>(); }, _sel);
}

template <typename N>
inline JoinGenomicCoords::iterator<N>::iterator(iterator::HicIt it,
                                                std::shared_ptr<const BinTable> bins)
    : _it(std::in_place_type_t<iterator::HicIt>(), std::move(it)), _bins(std::move(bins)) {
  std::ignore = **this;
}
template <typename N>
inline JoinGenomicCoords::iterator<N>::iterator(iterator::CoolerIt it,
                                                std::shared_ptr<const BinTable> bins)
    : _it(std::in_place_type_t<iterator::CoolerIt>(), std::move(it)), _bins(std::move(bins)) {
  std::ignore = **this;
}

template <typename N>
inline bool JoinGenomicCoords::iterator<N>::operator==(const iterator& other) const noexcept {
  return std::visit([&](const auto& it) {
    using T = decltype(it);
    assert(std::holds_alternative<T>(other));
    auto& other_ = std::get<T>(other);

    return it == other_;
  });
}
template <typename N>
inline bool JoinGenomicCoords::iterator<N>::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto JoinGenomicCoords::iterator<N>::operator*() const -> const_reference {
  assert(_bins);

  const auto& tp = std::visit([&](const auto& it) { return *it; });
  _value = Pixel<N>{_bins->at(tp.bin1_id), _bins->at(tp.bin2_id), tp.count};
  return _value;
}
template <typename N>
inline auto JoinGenomicCoords::iterator<N>::operator->() const -> const_pointer {
  return &(**this);
}

template <typename N>
inline auto JoinGenomicCoords::iterator<N>::operator++() -> iterator& {
  return std::visit([&](const auto& it) { return ++it; });
}
template <typename N>
inline auto JoinGenomicCoords::iterator<N>::operator++(int) -> iterator {
  return std::visit([&](const auto& it) { return it++; });
}
}  // namespace hictk::transformers
