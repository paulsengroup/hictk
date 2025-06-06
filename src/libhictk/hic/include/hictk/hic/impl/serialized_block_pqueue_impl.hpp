// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "hictk/common.hpp"

namespace hictk::hic::internal {
template <typename BlockID>
inline SerializedBlockPQueue<BlockID>::Record::operator bool() const noexcept {
  return status == Status::SUCCESS;
}

template <typename BlockID>
template <typename It>
inline SerializedBlockPQueue<BlockID>::SerializedBlockPQueue(It first_bid, It last_bid,
                                                             std::size_t producers_,
                                                             std::size_t capacity_)
    : _block_ids(first_bid, last_bid), _producers(producers_) {
  if (_block_ids.empty()) {
    _capacity = 0;
    return;
  }

  _capacity =  // NOLINTNEXTLINE(*-avoid-magic-numbers)
      std::clamp(capacity_ == 0 ? 3 * _producers : capacity_, std::size_t{2}, std::size_t{32});
  std::sort(_block_ids.begin(), _block_ids.end(), std::greater{});

  SPDLOG_DEBUG(FMT_STRING("initialized a BlockPQueue with capacity blocks with bids {}-{}"),
               _capacity, _block_ids.back(), _block_ids.front());
}

template <typename BlockID>
inline std::size_t SerializedBlockPQueue<BlockID>::size() const {
  [[maybe_unused]] const auto lck = std::scoped_lock(_mtx);
  return _buff.size();
}

template <typename BlockID>
inline std::size_t SerializedBlockPQueue<BlockID>::capacity() const noexcept {
  return _capacity;
}

template <typename BlockID>
inline bool SerializedBlockPQueue<BlockID>::try_enqueue(const BlockID& block_id,
                                                        const std::string& serialized_block,
                                                        std::chrono::milliseconds timeout) {
  const auto expiration = std::chrono::steady_clock::now() + timeout;
  while (HICTK_LIKELY(std::chrono::steady_clock::now() < expiration)) {
    {
      [[maybe_unused]] const auto lck = std::scoped_lock(_mtx);
      assert(!_block_ids.empty());
      assert(block_id == _block_ids.back() || block_id > _block_ids.back());
      if (HICTK_LIKELY(_buff.size() < _capacity - 1) || block_id == _block_ids.back()) {
        [[maybe_unused]] const auto inserted = _buff.emplace(block_id, serialized_block).second;
        assert(inserted);
        SPDLOG_DEBUG(
            FMT_STRING("SerializedBlockPQueue::try_enqueue(): successfully enqueued block {}"),
            block_id);
        return true;
      }
    }
    SPDLOG_DEBUG(FMT_STRING("SerializedBlockPQueue::try_enqueue(): failed to enqueue block {} "
                            "(queue is full). Sleeping "
                            "before trying one more time..."),
                 block_id);
    // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::this_thread::sleep_for(timeout / 25);
  }

  SPDLOG_DEBUG(
      FMT_STRING("SerializedBlockPQueue::try_enqueue(): failed to enqueue block {}: timed out!"),
      block_id);
  return false;
}

template <typename BlockID>
inline auto SerializedBlockPQueue<BlockID>::dequeue_timed(std::chrono::milliseconds timeout)
    -> Record {
  const auto expiration = std::chrono::steady_clock::now() + timeout;
  while (HICTK_LIKELY(std::chrono::steady_clock::now() < expiration)) {
    {
      [[maybe_unused]] const auto lck = std::scoped_lock(_mtx);
      auto result = dequeue_unsafe();
      if (result.status != Record::Status::NOT_AVAILABLE) {
        return result;
      }
    }
    SPDLOG_DEBUG(
        "SerializedBlockPQueue::dequeue_timed(): queue is empty. Sleeping before trying one more "
        "time...");
    // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::this_thread::sleep_for(timeout / 25);
  }

  SPDLOG_DEBUG("SerializedBlockPQueue::dequeue_timed(): operation timed out");
  return {{}, "", Record::Status::TIMEOUT};
}

template <typename BlockID>
inline void SerializedBlockPQueue<BlockID>::dequeue(std::vector<Record>& buffer) {
  buffer.clear();
  [[maybe_unused]] const auto lck = std::scoped_lock(_mtx);
  const auto fill_rate = static_cast<double>(_buff.size()) / static_cast<double>(_capacity);
  if (_block_ids.size() != _buff.size() && fill_rate < 0.5) {  // NOLINT(*-avoid-magic-numbers)
    // Queue is not full enough to be worth returning available items
    SPDLOG_DEBUG(FMT_STRING("SerializedBlockPQueue::dequeue(): not bothering dequeuing blocks "
                            "(queue is only {}% full)"),
                 100.0 * fill_rate);
    return;
  }
  dequeue_unsafe(buffer);
}

template <typename BlockID>
inline auto SerializedBlockPQueue<BlockID>::dequeue_unsafe() noexcept -> Record {
  if (HICTK_UNLIKELY(_block_ids.empty())) {
    SPDLOG_DEBUG(
        "SerializedBlockPQueue::dequeue_unsafe(): caught attempt to fetch block from a closed "
        "queue");
    return {{}, "", Record::Status::QUEUE_IS_CLOSED};
  }

  const auto wanted_bid = _block_ids.back();
  if (!_buff.empty()) {
    assert(!_block_ids.empty());
    auto& [bid, value] = *_buff.begin();
    assert(bid > wanted_bid || bid == wanted_bid);
    if (bid == wanted_bid) {
      assert(!value.empty());
      SPDLOG_DEBUG(FMT_STRING("SerializedBlockPQueue::dequeue_unsafe(): returning block {}..."),
                   bid);
      Record record{bid, std::move(value), Record::Status::SUCCESS};
      _buff.erase(_buff.begin());
      _block_ids.pop_back();
      return record;
    }
  }

  SPDLOG_DEBUG(
      FMT_STRING("SerializedBlockPQueue::dequeue_unsafe(): block {} has not yet been enqueued!"),
      wanted_bid);
  return {{}, "", Record::Status::NOT_AVAILABLE};
}

template <typename BlockID>
inline void SerializedBlockPQueue<BlockID>::dequeue_unsafe(std::vector<Record>& buffer) {
  buffer.clear();
  while (true) {
    auto record = dequeue_unsafe();
    if (!record) {
      SPDLOG_DEBUG(FMT_STRING("SerializedBlockPQueue::dequeue_unsafe(): dequeued {} blocks"),
                   buffer.size());
      return;
    }
    buffer.emplace_back(std::move(record));
  }
}

}  // namespace hictk::hic::internal
