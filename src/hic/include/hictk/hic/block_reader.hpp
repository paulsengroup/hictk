// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <limits>
#include <memory>

#include "hictk/chromosome.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/hic_file_stream.hpp"
#include "hictk/hic/index.hpp"

namespace hictk::hic::internal {

class BlockGrid {
  struct Node {
    std::shared_ptr<const BlockIndex> block{};
    const Node* next_right{};
    const Node* next_down{};
    std::size_t
        row{};  // These should be in absolute term (i.e. first/last row and col mapping to block)
    std::size_t
        col{};  // These should be in absolute term (i.e. first/last row and col mapping to block)
  };

  std::vector<Node> _grid{};

 public:
  class iterator;
  BlockGrid() = default;
  BlockGrid(const std::vector<BlockIndex>& index, std::size_t block_column_count);

  [[nodiscard]] auto begin() const noexcept -> iterator;
  [[nodiscard]] auto end() const noexcept -> iterator;
  [[nodiscard]] std::size_t size() const noexcept;

  class iterator {
    const Node* _node{};
    std::size_t _current_row{(std::numeric_limits<std::size_t>::max)()};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Node;
    using pointer = Node*;
    using const_pointer = const Node*;
    using reference = Node&;
    using const_reference = const Node&;
    using iterator_category = std::input_iterator_tag;

    iterator() = default;
    explicit iterator(const Node& head);
    [[nodiscard]] bool operator==(const iterator& other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator& other) const noexcept;

    auto operator++() noexcept -> iterator&;
    auto operator++(int) noexcept -> iterator;

    [[nodiscard]] auto operator*() noexcept -> const_reference;
    [[nodiscard]] auto operator->() noexcept -> const_pointer;
  };
};

class BinaryBuffer {
  std::string _buffer{};
  std::size_t _i{};

 public:
  BinaryBuffer() = default;
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type* = nullptr>
  T read();

  // Return the offset of the underlying buffer. Useful for error checking
  [[nodiscard]] std::size_t operator()() const noexcept;

  // Reset and return ref to underlying buffer so that buff can be refilled
  std::string& reset() noexcept;
};

class HiCBlockReader {
  std::shared_ptr<HiCFileStream> _hfs{};
  Index _index{};
  std::shared_ptr<BlockLRUCache> _blk_cache{};  // This should be passed in by file. Key should be
                                                // changed from size_t to {chrom1, chrom2, size_t}
  // We need the entire bin table in order to map pixels to abs bin ids
  std::shared_ptr<const BinTable> _bins{};
  BlockGrid _block_grid{};

  BinaryBuffer _bbuffer{};
  std::vector<SerializedPixel> _tmp_buffer{};

 public:
  HiCBlockReader() = default;
  HiCBlockReader(std::shared_ptr<HiCFileStream> hfs, const HiCFooter& footer,
                 std::shared_ptr<const BinTable> bins_, std::shared_ptr<BlockLRUCache> block_cache_,
                 const PixelCoordinates& coords1, const PixelCoordinates& coords2);

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] const BlockGrid& grid() const;

  [[nodiscard]] double sum() const noexcept;
  [[nodiscard]] double avg() const noexcept;

  [[nodiscard]] std::shared_ptr<const InteractionBlock> read(const BlockIndex& idx);

 private:
  void find_overlapping_blocks(const PixelCoordinates& coords1, const PixelCoordinates& coords2);
  [[nodiscard]] static Index read_index(HiCFileStream& hfs, const HiCFooter& footer);
  static void read_dispatcher_type1_block(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                          std::int32_t bin1Offset, std::int32_t bin2Offset,
                                          BinaryBuffer& src,
                                          std::vector<SerializedPixel>& dest) noexcept;
  template <typename Bin1Type, typename Bin2Type, typename CountType>
  static void read_type1_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<SerializedPixel>& dest) noexcept;

  template <typename CountType>
  static void read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<SerializedPixel>& dest) noexcept;
};

}  // namespace hictk::hic::internal

#include "../../../block_reader_impl.hpp"
