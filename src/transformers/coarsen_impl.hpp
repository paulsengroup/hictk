// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::transformers {

template <typename PixelIt>
inline CoarsenPixels<PixelIt>::CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel,
                                             std::size_t factor)
    : _first(std::move(first_pixel)), _last(std::move(last_pixel)), _factor(factor) {
  if (factor < 2) {
    throw std::logic_error("coarsening factor should be > 1");
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::begin() const -> iterator {
  return iterator{_first, _last, _factor};
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cbegin() const -> iterator {
  return this->begin();
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::end() const -> iterator {
  return iterator::at_end(_last, _factor);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cend() const -> iterator {
  return this->end();
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::read_all() const -> std::vector<ThinPixel<N>> {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<ThinPixel<N>> buff{};
  std::copy(this->begin(), this->end(), std::back_inserter(buff));
  return buff;
}

template <typename PixelIt>
[[nodiscard]] inline bool CoarsenPixels<PixelIt>::iterator::PixelCmp::operator()(
    const ThinPixel<N> &p1, const ThinPixel<N> &p2) const noexcept {
  if (p1.bin1_id != p2.bin1_id) {
    return p1.bin1_id < p2.bin1_id;
  }
  return p1.bin2_id < p2.bin2_id;
}

template <typename PixelIt>
inline CoarsenPixels<PixelIt>::iterator::iterator(PixelIt first, PixelIt last, std::size_t factor)
    : _pixel_it(std::move(first)),
      _pixel_last(std::move(last)),
      _merger(std::make_shared<PixelMerger>()),
      _factor(factor) {
  assert(_factor > 1);

  if (_pixel_it != _pixel_last) {
    _bin1_id_end = (_pixel_it->bin1_id / _factor) * _factor;
    process_next_row();
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::at_end(PixelIt last, std::size_t factor) -> iterator {
  iterator it{};
  it._pixel_it = last;
  it._pixel_last = last;
  it._factor = factor;
  return it;
}

template <typename PixelIt>
inline void CoarsenPixels<PixelIt>::iterator::process_next_row() {
  if (_merger.use_count() != 1) {
    _merger = std::make_shared<PixelMerger>();
  }
  _merger->clear();

  _bin1_id_end += _factor;

  while (_pixel_it != _pixel_last) {
    if (_pixel_it->bin1_id >= _bin1_id_end) {
      if (_merger->empty()) {
        process_next_row();  // found an empty row
      }
      break;
    }
    // Generate coarse pixel
    const auto pixel =
        ThinPixel<N>{_pixel_it->bin1_id / _factor, _pixel_it->bin2_id / _factor, _pixel_it->count};

    // Merge or insert coarse pixel
    auto it = _merger->find(pixel);
    if (it != _merger->end()) {
      it->count += pixel.count;
    } else {
      _merger->emplace(std::move(pixel));
    }
    ++_pixel_it;
  }

  _it = _merger->begin();
  if (_merger->empty()) {
    assert(_pixel_it == _pixel_last);
    *this = at_end(_pixel_last, _factor);
  }
}

template <typename PixelIt>
inline bool CoarsenPixels<PixelIt>::iterator::operator==(
    const CoarsenPixels::iterator &other) const noexcept {
  return _merger == other._merger && _it == other._it && this->_factor == other._factor;
}

template <typename PixelIt>
inline bool CoarsenPixels<PixelIt>::iterator::operator!=(
    const CoarsenPixels::iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator*() const -> const_reference {
  assert(this->_merger);
  assert(this->_it != this->_merger->end());
  return *this->_it;
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator->() const -> const_pointer {
  return &(**this);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator++() -> iterator & {
  assert(_merger);

  if (++this->_it == this->_merger->end()) {
    this->process_next_row();
  }
  return *this;
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace hictk::transformers
