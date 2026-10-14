// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#if __has_include(<readerwriterqueue.h>)
#include <readerwriterqueue.h>
#else
#include <readerwriterqueue/readerwriterqueue.h>
#endif

#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {  // NOLINT(*-concat-nested-namespaces)

template <typename N>
void MatrixInteractionBlockFlat<N>::emplace_back(Pixel<N> p) {
  emplace_back(p.to_thin());
}

template <typename N>
void MatrixInteractionBlockFlat<N>::emplace_back(const ThinPixel<N> &p) {
  bin1_ids.push_back(p.bin1_id);
  bin2_ids.push_back(p.bin2_id);
  counts.push_back(p.count);
}

template <typename N>
std::size_t MatrixInteractionBlockFlat<N>::size() const noexcept {
  return bin1_ids.size();
}

template <typename N>
std::string MatrixInteractionBlockFlat<N>::serialize(BinaryBuffer &buffer, ZSTD_CCtx_s &compressor,
                                                     std::string &compression_buffer,
                                                     int compression_lvl, bool clear) const {
  if (size() == 0) {
    return "";
  }

  if (clear) {
    buffer.clear();
  }

  buffer.write(bin1_ids);
  buffer.write(bin2_ids);
  buffer.write(counts);

  const auto buff_size = ZSTD_compressBound(buffer.get().size() * sizeof(char));
  compression_buffer.resize(buff_size);

  const std::size_t compressed_size = ZSTD_compressCCtx(
      &compressor, static_cast<void *>(compression_buffer.data()),
      compression_buffer.size() * sizeof(char), static_cast<const void *>(buffer.get().data()),
      buffer.get().size() * sizeof(char), compression_lvl);
  if (ZSTD_isError(compressed_size) != 0) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  compression_buffer.resize(compressed_size);

  buffer.clear();
  buffer.write(size());
  buffer.write(compression_buffer, false);

  return buffer.get();
}

template <typename N>
[[nodiscard]] std::vector<ThinPixel<N>> MatrixInteractionBlockFlat<N>::deserialize(
    BinaryBuffer &buffer, ZSTD_DCtx_s &decompressor, std::string &decompression_buffer) {
  const auto size_ = buffer.read<std::size_t>();
  std::vector<ThinPixel<N>> pixels(size_);

  const auto decompressed_size =
      size_ * (sizeof(std::uint64_t) + sizeof(std::uint64_t) + sizeof(N));
  decompression_buffer.resize(decompressed_size);

  const auto compressed_buffer = std::string_view{buffer.get()}.substr(sizeof(size_));

  const auto status = ZSTD_decompressDCtx(
      &decompressor, decompression_buffer.data(), decompression_buffer.size() * sizeof(char),
      compressed_buffer.data(), compressed_buffer.size() * sizeof(char));

  if (ZSTD_isError(status) != 0) {
    throw std::runtime_error(ZSTD_getErrorName(status));
  }
  buffer.clear();
  buffer.write(decompression_buffer);

  for (auto &p : pixels) {
    p.bin1_id = buffer.read<std::uint64_t>();
  }
  for (auto &p : pixels) {
    p.bin2_id = buffer.read<std::uint64_t>();
  }
  for (auto &p : pixels) {
    p.count = buffer.read<N>();
  }

  return pixels;
}

template <typename PixelIt, typename>
void HiCInteractionToBlockMapper::append_pixels(PixelIt first_pixel, const PixelIt &last_pixel,
                                                bool validate, std::uint32_t update_frequency) {
  using PixelT = remove_cvref_t<decltype(*first_pixel)>;
  static_assert(std::is_same_v<PixelT, ThinPixel<float>> || std::is_same_v<PixelT, Pixel<float>>);

  SPDLOG_DEBUG(FMT_STRING("mapping pixels to interaction blocks at resolution {}..."),
               _bin_table->resolution());

  [[maybe_unused]] auto t0 = std::chrono::steady_clock::now();
  for (std::size_t i = 0; first_pixel != last_pixel; ++i) {
    if (_pending_pixels >= _chunk_size) {
      write_blocks();
    }

    add_pixel(*first_pixel, validate);
    std::ignore = ++first_pixel;  // NOLINT(*-pointer-arithmetic)

    if constexpr (SPDLOG_ACTIVE_LEVEL >= SPDLOG_LEVEL_INFO) {
      if (i == update_frequency) {
        const auto t1 = std::chrono::steady_clock::now();
        const auto delta =
            static_cast<double>(
                std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
            1000.0;
        SPDLOG_INFO(FMT_STRING("ingesting pixels at {:.0f} pixels/s..."),
                    static_cast<double>(update_frequency) / delta);
        t0 = t1;
        i = 0;
      }
    }
  }
}

namespace internal {
// I am declaring this as a freestanding function instead that as a lambda to work around a compiler
// bug in MSVC
template <typename PixelT>
Pixel<float> process_pixel_interaction_block(const BinTable &bin_table, PixelT pixel) {
  static_assert(std::is_same_v<PixelT, ThinPixel<float>> || std::is_same_v<PixelT, Pixel<float>>);
  constexpr bool is_thin_pixel = std::is_same_v<PixelT, ThinPixel<float>>;

  if constexpr (is_thin_pixel) {
    return Pixel(bin_table, pixel);
  } else {
    return pixel;
  }
}
}  // namespace internal

template <typename PixelIt, typename>
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void HiCInteractionToBlockMapper::append_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                                bool validate, BS::light_thread_pool &tpool,
                                                std::uint32_t update_frequency) {
  if (tpool.get_thread_count() < 2) {
    return append_pixels(first_pixel, last_pixel, validate, update_frequency);
  }

  SPDLOG_DEBUG("mapping pixels to interaction blocks using 2 threads...");

  std::atomic<bool> early_return = false;
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  moodycamel::BlockingReaderWriterQueue<Pixel<float>> queue(10'000);

  auto writer = tpool.submit_task([&]() {
    try {
      [[maybe_unused]] auto t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; first_pixel != last_pixel && !early_return; ++i) {
        const auto pixel = internal::process_pixel_interaction_block(*_bin_table, *first_pixel);
        assert(pixel.count != 0);
        std::ignore = ++first_pixel;  // NOLINT(*-pointer-arithmetic)

        while (!queue.try_enqueue(pixel)) {
          if (early_return) {
            return;
          }
        }
        if constexpr (SPDLOG_ACTIVE_LEVEL >= SPDLOG_LEVEL_INFO) {
          if (i == update_frequency) {
            const auto t1 = std::chrono::steady_clock::now();
            const auto delta =
                static_cast<double>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
                1000.0;
            SPDLOG_INFO(FMT_STRING("ingesting pixels at {:.0f} pixels/s..."),
                        static_cast<double>(update_frequency) / delta);
            t0 = t1;
            i = 0;
          }
        }
      }
      queue.enqueue(Pixel<float>{});

    } catch (...) {
      early_return = true;
      throw;
    }
  });

  auto reader = tpool.submit_task([&]() {
    try {
      Pixel<float> p{};
      while (!early_return) {
        if (_pending_pixels >= _chunk_size) {
          write_blocks();
        }

        queue.wait_dequeue(p);
        if (p.coords.bin1.id() == Bin::null_id && p.coords.bin2.id() == Bin::null_id &&
            p.count == 0) {
          return;
        }

        add_pixel(p, validate);
        ++_processed_pixels;
        ++_pending_pixels;
      }
    } catch (...) {
      early_return = true;
      throw;
    }
  });

  writer.get();
  reader.get();
}

template <typename N>
auto HiCInteractionToBlockMapper::map(const ThinPixel<N> &p) const -> BlockID {
  return map(Pixel<N>(*_bin_table, p));
}

template <typename N>
auto HiCInteractionToBlockMapper::map(const Pixel<N> &p) const -> BlockID {
  const auto &bin1 = p.coords.bin1;
  const auto &bin2 = p.coords.bin2;

  const auto &chrom1 = bin1.chrom();
  const auto &chrom2 = bin2.chrom();

  const auto bin1_id = bin1.rel_id();
  const auto bin2_id = bin2.rel_id();

  const auto block_id = p.coords.is_intra()
                            ? _mappers_intra.at(chrom1)(bin1_id, bin2_id)
                            : _mappers_inter.at(std::make_pair(chrom1, chrom2))(bin1_id, bin2_id);

  return {chrom1.id(), chrom2.id(), block_id};
}

template <typename N>
void HiCInteractionToBlockMapper::add_pixel(const ThinPixel<N> &p, bool validate) {
  if (validate) {
    if (HICTK_UNLIKELY(p.count == 0)) {
      throw std::runtime_error(fmt::format(FMT_STRING("({}) found a pixel of value 0"), p));
    }
    if (HICTK_UNLIKELY(p.bin1_id > _bin_table->size())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"), p, p.bin1_id));
    }
    if (HICTK_UNLIKELY(p.bin2_id > _bin_table->size())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"), p, p.bin2_id));
    }
    if (HICTK_UNLIKELY(p.bin1_id > p.bin2_id)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("({}) bin1_id is greater than bin2_id: {} > {}"), p, p.bin1_id, p.bin2_id));
    }
  }

  add_pixel(Pixel(*_bin_table, p), false);
}

template <typename N>
void HiCInteractionToBlockMapper::add_pixel(const Pixel<N> &p, bool validate) {
  if (validate) {
    if (HICTK_UNLIKELY(p.count == 0)) {
      throw std::runtime_error(fmt::format(FMT_STRING("({}) found a pixel of value 0"), p));
    }

    if (HICTK_UNLIKELY(!chromosomes().contains(p.coords.bin1.chrom().id()))) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("({}) invalid chromosome id {}"), p, p.coords.bin1.chrom().id()));
    }

    if (HICTK_UNLIKELY(!chromosomes().contains(p.coords.bin2.chrom().id()))) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("({}) invalid chromosome id {}"), p, p.coords.bin2.chrom().id()));
    }

    if (const auto bin_id = p.coords.bin1.id(); HICTK_UNLIKELY(bin_id > _bin_table->size())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"), p, bin_id));
    }

    if (const auto bin_id = p.coords.bin2.id(); HICTK_UNLIKELY(bin_id > _bin_table->size())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("({}) invalid bin id {}: bin maps outside of the bin table"), p, bin_id));
    }

    if (HICTK_UNLIKELY(p.coords.bin1.id() > p.coords.bin2.id())) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("({}) bin1_id is greater than bin2_id: {} > {}"), p,
                      p.coords.bin1.id(), p.coords.bin2.id()));
    }
  }

  const auto bid = map(p);

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();
  const auto chrom_pair = std::make_pair(chrom1, chrom2);

  auto match1 = _blocks.find(bid);
  if (match1 != _blocks.end()) {
    _pixel_sums.at(chrom_pair) += p.count;
    match1->second.emplace_back(p.to_thin());
  } else {
    _pixel_sums.emplace(chrom_pair, p.count);
    auto [it, _] = _blocks.emplace(bid, MatrixInteractionBlockFlat<float>{});
    it->second.emplace_back(p.to_thin());
  }

  auto match2 = _chromosome_index.find(chrom_pair);
  if (match2 != _chromosome_index.end()) {
    match2->second.emplace(bid);
  } else {
    _chromosome_index.emplace(chrom_pair, phmap::btree_set<BlockID>{bid});
  }
  ++_processed_pixels;
  ++_pending_pixels;
}

}  // namespace hictk::hic::internal
