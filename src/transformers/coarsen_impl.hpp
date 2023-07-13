// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::transformers {

template <typename PixelIt>
inline CoarsenPixels<PixelIt>::CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel,
                                             std::shared_ptr<const BinTable> source_bins,
                                             std::size_t factor)
    : _first(std::move(first_pixel)),
      _last(std::move(last_pixel)),
      _src_bins(std::move(source_bins)),
      _dest_bins(std::make_shared<const BinTable>(_src_bins->chromosomes(),
                                                  _src_bins->bin_size() * factor)),
      _factor(factor) {
  if (factor < 2) {
    throw std::logic_error("coarsening factor should be > 1");
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::begin() const -> iterator {
  return iterator{_first, _last, _src_bins, _dest_bins};
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cbegin() const -> iterator {
  return this->begin();
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::end() const -> iterator {
  return iterator::at_end(_last, _src_bins, _dest_bins);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cend() const -> iterator {
  return this->end();
}

template <typename PixelIt>
inline const BinTable &CoarsenPixels<PixelIt>::src_bins() const noexcept {
  return *this->src_bins_ptr();
}
template <typename PixelIt>
inline const BinTable &CoarsenPixels<PixelIt>::dest_bins() const noexcept {
  return *this->dest_bins_ptr();
}
template <typename PixelIt>
inline std::shared_ptr<const BinTable> CoarsenPixels<PixelIt>::src_bins_ptr() const noexcept {
  return this->_src_bins;
}
template <typename PixelIt>
inline std::shared_ptr<const BinTable> CoarsenPixels<PixelIt>::dest_bins_ptr() const noexcept {
  return this->_dest_bins;
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
inline CoarsenPixels<PixelIt>::iterator::iterator(PixelIt first, PixelIt last,
                                                  std::shared_ptr<const BinTable> src_bins,
                                                  std::shared_ptr<const BinTable> dest_bins)
    : _pixel_it(std::move(first)),
      _pixel_last(std::move(last)),
      _src_bins(std::move(src_bins)),
      _dest_bins(std::move(dest_bins)),
      _merger(std::make_shared<PixelMerger>()) {
  assert(_dest_bins->bin_size() > _src_bins->bin_size());

  if (_pixel_it == _pixel_last) {
    *this = at_end(_pixel_last, _src_bins, _dest_bins);
  } else {
    process_next_row();
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::at_end(PixelIt last,
                                                     std::shared_ptr<const BinTable> src_bins,
                                                     std::shared_ptr<const BinTable> dest_bins)
    -> iterator {
  iterator it{};
  it._pixel_it = last;
  it._pixel_last = last;
  it._src_bins = std::move(src_bins);
  it._dest_bins = std::move(dest_bins);
  return it;
}

template <typename PixelIt>
inline void CoarsenPixels<PixelIt>::iterator::process_next_row() {
  if (_merger.use_count() != 1) {
    _merger = std::make_shared<PixelMerger>();
  }
  _merger->clear();

  const auto factor = _dest_bins->bin_size() / _src_bins->bin_size();
  _bin1_id_end += factor;

  while (_pixel_it != _pixel_last) {
    // We need to map pixel coordinates instead of just dividing bin ids by the coarsening factor to
    // avoid mapping the last pixel in a chromosome i and the first pixel in chromosome i+1 as to
    // the same coarse bin
    const PixelCoordinates src_coords{_src_bins->at(_pixel_it->bin1_id),
                                      _src_bins->at(_pixel_it->bin2_id)};
    const PixelCoordinates dest_coords{_dest_bins->at(src_coords.bin1.interval()).first,
                                       _dest_bins->at(src_coords.bin2.interval()).first};

    if (src_coords.bin1.rel_id() >= _bin1_id_end) {
      if (_merger->empty()) {
        process_next_row();  // found an empty row
      }
      break;
    }

    auto pixel = ThinPixel<N>{dest_coords.bin1.id(), dest_coords.bin2.id(), _pixel_it->count};

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
    *this = at_end(_pixel_last, _src_bins, _dest_bins);
  }
}

template <typename PixelIt>
inline bool CoarsenPixels<PixelIt>::iterator::operator==(
    const CoarsenPixels::iterator &other) const noexcept {
  return _merger == other._merger && _it == other._it && this->_src_bins == other._src_bins &&
         this->_dest_bins == other._dest_bins;
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
