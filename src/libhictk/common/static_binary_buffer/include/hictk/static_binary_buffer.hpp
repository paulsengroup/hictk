// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <type_traits>

namespace hictk::internal {

template <typename... Ts>
class StaticBinaryBuffer {
  static_assert((std::is_arithmetic_v<Ts> && ...), "Not all types are arithmetic!");

  // We are allocating one extra byte so that we can store '\0' as the last value
  // (just in case the string_view returned by operator()() is passed to some function
  // that expects null-terminated strings)
  static constexpr std::size_t SIZE = (sizeof(Ts) + ... + 1);

  alignas(std::uint64_t) std::array<char, SIZE> _buff{'\0'};

 public:
  StaticBinaryBuffer() = default;
  explicit StaticBinaryBuffer(const Ts &...inputs) noexcept {
    std::size_t i = 0;
    (
        [&] {
          // This is just a best-effort attempt to make the compiler do the data copying at compile
          // time.
          // When this fails, and the copying is done at run time, this implementation will be much
          // slower than e.g. std::memcopy.
          // However, StaticBinaryBuffer should never be used to more than say a few dozens of bytes
          // worth of data, so even if the data copying ends up being slow, it shouldn't be a big
          // deal.
          // NOLINTNEXTLINE(*-type-reinterpret-cast)
          const auto *src = reinterpret_cast<const char *>(&inputs);
          for (std::size_t j = 0; j < sizeof(inputs); ++j, ++i) {
            _buff[i] = *(src + j);  // NOLINT(*-bounds-pointer-arithmetic)
          }
        }(),
        ...);

    assert(!_buff.empty());
    _buff.back() = '\0';
  }

  constexpr std::string_view operator()() const noexcept {
    return std::string_view{_buff.data(), _buff.size() - 1};
  }
};

}  // namespace hictk::internal
