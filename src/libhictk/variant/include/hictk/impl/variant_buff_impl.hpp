// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/generic_variant.hpp"

namespace hictk::internal {

// NOLINTBEGIN(modernize-use-nodiscard)

template <typename T>
inline VariantBuffer::VariantBuffer(std::vector<T> data) : _buff(std::move(data)) {}
template <typename InputIt>
inline VariantBuffer::VariantBuffer(InputIt first, InputIt last)
    : VariantBuffer(std::vector<typename InputIt::value_type>(first, last)) {}
template <typename T>
inline VariantBuffer::VariantBuffer(std::size_t size, T default_value)
    : VariantBuffer(std::vector<T>(size, default_value)) {}

template <typename T>
inline VariantBuffer &VariantBuffer::operator=(std::vector<T> buff) noexcept {
  _buff = std::move(buff);

  return *this;
}

template <typename T>
inline typename std::vector<T>::iterator VariantBuffer::VariantBuffer::begin() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).begin();
}
template <typename T>
inline typename std::vector<T>::iterator VariantBuffer::VariantBuffer::end() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).end();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::begin() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).begin();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::end() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).end();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::cbegin() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).begin();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::cend() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).end();
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline std::size_t VariantBuffer::VariantBuffer::size() const noexcept {
  assert(!_buff.valueless_by_exception());
  std::size_t size{};
  std::visit([&](const auto &buff) { size = buff.size(); }, _buff);
  return size;
}

template <typename T>
inline std::size_t VariantBuffer::VariantBuffer::size() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).size();
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline std::size_t VariantBuffer::VariantBuffer::capacity() const noexcept {
  assert(!_buff.valueless_by_exception());
  std::size_t capacity{};
  std::visit([&](const auto &buff) { capacity = buff.capacity(); }, _buff);
  return capacity;
}

template <typename T>
inline std::size_t VariantBuffer::VariantBuffer::capacity() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).capacity();
}

template <typename T>
inline void VariantBuffer::VariantBuffer::reserve(std::size_t new_size) {
  assert(!_buff.valueless_by_exception());
  std::get<std::vector<T>>(_buff).reserve(new_size);
}

template <typename T>
inline void VariantBuffer::VariantBuffer::resize(std::size_t new_size) {
  assert(!_buff.valueless_by_exception());
  std::get<std::vector<T>>(_buff).resize(new_size);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline bool VariantBuffer::VariantBuffer::empty() const noexcept { return size() == 0; }

template <typename T>
inline bool VariantBuffer::VariantBuffer::empty() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).empty();
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline void VariantBuffer::VariantBuffer::clear() noexcept {
  assert(!_buff.valueless_by_exception());
  std::visit([](auto &buff) { buff.clear(); }, _buff);
}

template <typename T>
inline void VariantBuffer::VariantBuffer::clear() {
  assert(!_buff.valueless_by_exception());
  std::get<std::vector<T>>(_buff).clear();
}

template <typename T>
[[nodiscard]] inline T &VariantBuffer::at(std::size_t i) {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).at(i);
}
template <typename T>
[[nodiscard]] inline const T &VariantBuffer::at(std::size_t i) const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).at(i);
}

[[nodiscard]] inline GenericVariant VariantBuffer::at(std::size_t i) const {
  assert(!_buff.valueless_by_exception());
  GenericVariant n{};
  std::visit([&](const auto &buff) { n = buff.at(i); }, _buff);
  return n;
}

[[nodiscard]] inline GenericVariant VariantBuffer::operator[](std::size_t i) const {
  assert(!_buff.valueless_by_exception());
  GenericVariant n{};
  std::visit(
      [&](const auto &buff) {
        assert(i < buff.size());
        n = buff[i];
      },
      _buff);
  return n;
}

template <typename T>
inline T &VariantBuffer::front() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).front();
}
template <typename T>
inline const T &VariantBuffer::front() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).front();
}

template <typename T>
inline T &VariantBuffer::back() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).back();
}
template <typename T>
inline const T &VariantBuffer::back() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).back();
}

template <typename T>
inline T *VariantBuffer::data() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).data();
}
template <typename T>
inline const T *VariantBuffer::data() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff).data();
}

template <typename T>
inline std::vector<T> &VariantBuffer::get() {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff);
}
template <typename T>
inline const std::vector<T> &VariantBuffer::get() const {
  assert(!_buff.valueless_by_exception());
  return std::get<std::vector<T>>(_buff);
}

constexpr auto VariantBuffer::get() -> BuffT & { return _buff; }
constexpr auto VariantBuffer::get() const -> const BuffT & { return _buff; }

template <typename T>
bool VariantBuffer::holds_alternative() const noexcept {
  assert(!_buff.valueless_by_exception());
  return std::holds_alternative<std::vector<T>>(_buff);
}

// NOLINTEND(modernize-use-nodiscard)

}  // namespace hictk::internal
