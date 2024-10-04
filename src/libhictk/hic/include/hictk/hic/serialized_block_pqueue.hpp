// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

namespace hictk::hic::internal {

template <typename BlockID>
class SerializedBlockPQueue {
  std::vector<BlockID> _block_ids{};
  phmap::btree_map<BlockID, std::string> _buff{};
  mutable std::mutex _mtx{};
  std::size_t _capacity{256};
  std::size_t _producers{1};

 public:
  struct Record {
    enum class Status : std::uint_fast8_t { SUCCESS, TIMEOUT, NOT_AVAILABLE, QUEUE_IS_CLOSED };

    BlockID bid{};
    std::string serialized_block{};
    Status status{Status::SUCCESS};

    [[nodiscard]] explicit operator bool() const noexcept;
  };

  SerializedBlockPQueue() = default;
  template <typename It>
  SerializedBlockPQueue(It first_bid, It last_bid, std::size_t producers,
                        std::size_t capacity_ = 0);

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] std::size_t capacity() const noexcept;

  [[nodiscard]] bool try_enqueue(const BlockID& block_id, const std::string& serialized_block,
                                 std::chrono::milliseconds timeout = std::chrono::seconds(1));
  [[nodiscard]] auto dequeue_timed(
      std::chrono::milliseconds timeout = std::chrono::milliseconds(50)) -> Record;
  void dequeue(std::vector<Record>& buffer);

 private:
  [[nodiscard]] auto dequeue_unsafe() noexcept -> Record;
};

}  // namespace hictk::hic::internal

#include "./impl/serialized_block_pqueue_impl.hpp"
