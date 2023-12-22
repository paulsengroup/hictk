// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator<(const Node &other) const noexcept {
  assert(!!pixel);
  assert(!!other.pixel);
  if (pixel.bin1_id != other.pixel.bin1_id) {
    return pixel.bin1_id < other.pixel.bin1_id;
  }
  return pixel.bin2_id < other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator>(const Node &other) const noexcept {
  assert(!!pixel);
  assert(!!other.pixel);
  if (pixel.bin1_id != other.pixel.bin1_id) {
    return pixel.bin1_id > other.pixel.bin1_id;
  }
  return pixel.bin2_id > other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator==(const Node &other) const noexcept {
  return pixel.bin1_id == other.pixel.bin1_id && pixel.bin2_id == other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator!=(const Node &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline PixelMerger<PixelIt>::PixelMerger(std::vector<PixelIt> heads, std::vector<PixelIt> tails)
    : PixelMerger(heads.begin(), heads.end(), tails.begin()) {
  assert(heads.size() == tails.size());
}

template <typename PixelIt>
template <typename ItOfPixelIt>
inline PixelMerger<PixelIt>::PixelMerger(ItOfPixelIt first_head, ItOfPixelIt last_head,
                                         ItOfPixelIt first_tail) {
  while (first_head != last_head) {
    if (*first_head != *first_tail) {
      _heads.emplace_back(std::move(*first_head));
      _tails.emplace_back(std::move(*first_tail));
    }
    ++first_head;
    ++first_tail;
  }
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::begin() const -> iterator {
  return iterator{_heads, _tails};
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::end() const noexcept -> iterator {
  return {};
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::read_all() const -> std::vector<PixelT> {
  std::vector<PixelT> buff{};
  std::copy(begin(), end(), std::back_inserter(buff));
  return buff;
}

template <typename PixelIt>
inline PixelMerger<PixelIt>::iterator::iterator(const std::vector<PixelIt> &heads,
                                                const std::vector<PixelIt> &tails)
    : _pqueue(std::make_shared<PQueueT>()),
      _heads(std::make_shared<std::vector<PixelIt>>(heads)),
      _tails(std::make_shared<std::vector<PixelIt>>(tails)) {
  assert(heads.size() == tails.size());
  for (auto &it : *_heads) {
    _pqueue->emplace(Node{*it++, _pqueue->size()});
  }
  _value = next();
}

template <typename PixelIt>
inline PixelMerger<PixelIt>::iterator::iterator(const iterator &other)
    : _value(other._value),
      _pqueue(other._pqueue ? std::make_shared<PQueueT>(*other._pqueue) : nullptr),
      _heads(other._heads ? std::make_shared<std::vector<PixelIt>>(*other._heads) : nullptr),
      _tails(other._tails ? std::make_shared<std::vector<PixelIt>>(*other._tails) : nullptr),
      _i(other._i) {}

template <typename PixelIt>
inline PixelMerger<PixelIt>::iterator::iterator(iterator &&other) noexcept
    : _value(other._value),
      _pqueue(std::move(other._pqueue)),
      _heads(std::move(other._heads)),
      _tails(std::move(other._tails)),
      _i(other._i) {}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::operator=(const iterator &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

  _value = other._value;
  _pqueue = other._pqueue ? std::make_shared<PQueueT>(*other._pqueue) : nullptr;
  _heads = other._heads ? std::make_shared<std::vector<PixelIt>>(*other._heads) : nullptr;
  _tails = other._tails ? std::make_shared<std::vector<PixelIt>>(*other._tails) : nullptr;
  _i = other._i;

  return *this;
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::operator=(iterator &&other) noexcept -> iterator & {
  if (this == &other) {
    return *this;
  }

  _value = std::move(other._value);
  _pqueue = std::move(other._pqueue);
  _heads = std::move(other._heads);
  _tails = std::move(other._tails);
  _i = other._i;

  return *this;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::iterator::operator==(const iterator &other) const noexcept {
  if (!_heads || !other._heads) {
    // check if we are at end
    return _heads == other._heads;
  }
  return _heads == other._heads && _i == other._i;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::operator*() const noexcept -> const ThinPixel<N> & {
  return _value;
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::operator->() const noexcept -> const ThinPixel<N> * {
  return &(**this);
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::operator++() -> iterator & {
  assert(_heads);
  _value = next();
  if (!_value) {
    _heads = nullptr;
    _tails = nullptr;
  }

  return *this;
}

template <typename PixelIt>
inline void PixelMerger<PixelIt>::iterator::replace_top_node() {
  assert(_heads);
  assert(_tails);
  const auto i = _pqueue->top().i;
  _pqueue->pop();
  if (auto &it = (*_heads)[i]; it != (*_tails)[i]) {
    _pqueue->emplace(Node{*it, i});
    ++it;
  }
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::iterator::next() -> PixelT {
  if (_pqueue->empty()) {
    return {};
  }

  auto current_node = _pqueue->top();
  replace_top_node();

  while (!_pqueue->empty()) {
    const auto next_node = _pqueue->top();
    if (next_node != current_node) {
      break;
    }
    current_node.pixel.count += next_node.pixel.count;
    replace_top_node();
  }
  ++_i;
  return current_node.pixel;
}

}  // namespace hictk::transformers
