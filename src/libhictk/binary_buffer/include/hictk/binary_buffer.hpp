// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstring>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace hictk {

class BinaryBuffer {
  std::string _buffer{};
  std::size_t _i{};

 public:
  BinaryBuffer() = default;

  // NOLINTNEXTLINE
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
  T read();
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
  void read(T& buff);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
  void read(std::vector<T>& buff);
  void read(std::string& buff, std::size_t n);
  void read(char* buff, std::size_t n);
  std::string getline(char delim = '\n');
  // NOLINTNEXTLINE
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
  void write(T data);
  void write(const std::string& data, bool add_nullterm = true);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>>* = nullptr>
  void write(const std::vector<T>& data);

  // Return the offset of the underlying buffer. Useful for error checking
  [[nodiscard]] std::size_t operator()() const noexcept;

  // Reset and return ref to underlying buffer so that buff can be refilled
  std::string& reset() noexcept;

  void clear() noexcept;

  [[nodiscard]] const std::string& get() const noexcept;
};

}  // namespace hictk

#include "./impl/binary_buffer_impl.hpp"  // NOLINT
