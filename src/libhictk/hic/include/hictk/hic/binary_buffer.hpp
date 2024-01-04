// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string>
#include <type_traits>

namespace hictk::hic::internal {

class BinaryBuffer {
  std::string _buffer{};
  std::size_t _i{};

 public:
  BinaryBuffer() = default;

  // NOLINTNEXTLINE
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type* = nullptr>
  T read();
  // NOLINTNEXTLINE
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type* = nullptr>
  void write(T data);
  void write(const std::string& data);

  // Return the offset of the underlying buffer. Useful for error checking
  [[nodiscard]] std::size_t operator()() const noexcept;

  // Reset and return ref to underlying buffer so that buff can be refilled
  std::string& reset() noexcept;

  void clear() noexcept;

  [[nodiscard]] const std::string& get() const noexcept;
};

}  // namespace hictk::hic::internal

#include "./impl/binary_buffer_impl.hpp"  // NOLINT
