// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <utility>
#include <variant>
#include <vector>

namespace coolerpp::internal {
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
  this->_buff = std::move(buff);

  return *this;
}

template <typename T>
inline typename std::vector<T>::iterator VariantBuffer::VariantBuffer::begin() {
  return std::get<std::vector<T>>(this->_buff).begin();
}
template <typename T>
inline typename std::vector<T>::iterator VariantBuffer::VariantBuffer::end() {
  return std::get<std::vector<T>>(this->_buff).end();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::begin() const {
  return std::get<std::vector<T>>(this->_buff).begin();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::end() const {
  return std::get<std::vector<T>>(this->_buff).end();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::cbegin() const {
  return std::get<std::vector<T>>(this->_buff).begin();
}
template <typename T>
inline typename std::vector<T>::const_iterator VariantBuffer::VariantBuffer::cend() const {
  return std::get<std::vector<T>>(this->_buff).end();
}

inline std::size_t VariantBuffer::VariantBuffer::size() const noexcept {
  std::size_t size{};
  std::visit([&](const auto &buff) { size = buff.size(); }, this->_buff);
  return size;
}

template <typename T>
inline std::size_t VariantBuffer::VariantBuffer::size() const {
  return std::get<std::vector<T>>(this->_buff).size();
}

inline std::size_t VariantBuffer::VariantBuffer::capacity() const noexcept {
  std::size_t capacity{};
  std::visit([&](const auto &buff) { capacity = buff.capacity(); }, this->_buff);
  return capacity;
}

template <typename T>
inline std::size_t VariantBuffer::VariantBuffer::capacity() const {
  return std::get<std::vector<T>>(this->_buff).capacity();
}

template <typename T>
inline void VariantBuffer::VariantBuffer::reserve(std::size_t new_size) {
  std::get<std::vector<T>>(this->_buff).reserve(new_size);
}

template <typename T>
inline void VariantBuffer::VariantBuffer::resize(std::size_t new_size) {
  std::get<std::vector<T>>(this->_buff).resize(new_size);
}

inline bool VariantBuffer::VariantBuffer::empty() const noexcept { return this->size() == 0; }

template <typename T>
inline bool VariantBuffer::VariantBuffer::empty() const {
  return std::get<std::vector<T>>(this->_buff).empty();
}

inline void VariantBuffer::VariantBuffer::clear() noexcept {
  std::visit([](auto &buff) { buff.clear(); }, this->_buff);
}

template <typename T>
inline void VariantBuffer::VariantBuffer::clear() noexcept {
  std::get<std::vector<T>>(this->_buff).clear();
}

template <typename T>
[[nodiscard]] inline T &VariantBuffer::at(std::size_t i) {
  return std::get<std::vector<T>>(this->_buff).at(i);
}
template <typename T>
[[nodiscard]] inline const T &VariantBuffer::at(std::size_t i) const {
  return std::get<std::vector<T>>(this->_buff).at(i);
}

[[nodiscard]] inline GenericVariant VariantBuffer::at(std::size_t i) const {
  GenericVariant n{};
  std::visit([&](const auto &buff) { n = buff.at(i); }, this->_buff);
  return n;
}

[[nodiscard]] inline GenericVariant VariantBuffer::operator[](std::size_t i) const {
  GenericVariant n{};
  std::visit(
      [&](const auto &buff) {
        assert(i < buff.size());
        n = buff[i];
      },
      this->_buff);
  return n;
}

template <typename T>
inline T &VariantBuffer::front() {
  return std::get<std::vector<T>>(this->_buff).front();
}
template <typename T>
inline const T &VariantBuffer::front() const {
  return std::get<std::vector<T>>(this->_buff).front();
}

template <typename T>
inline T &VariantBuffer::back() {
  return std::get<std::vector<T>>(this->_buff).back();
}
template <typename T>
inline const T &VariantBuffer::back() const {
  return std::get<std::vector<T>>(this->_buff).back();
}

template <typename T>
inline T *VariantBuffer::data() {
  return std::get<std::vector<T>>(this->_buff).data();
}
template <typename T>
inline const T *VariantBuffer::data() const {
  return std::get<std::vector<T>>(this->_buff).data();
}

template <typename T>
inline std::vector<T> &VariantBuffer::get() {
  return std::get<std::vector<T>>(this->_buff);
}
template <typename T>
inline const std::vector<T> &VariantBuffer::get() const {
  return std::get<std::vector<T>>(this->_buff);
}

constexpr auto VariantBuffer::get() -> BuffT & { return this->_buff; }
constexpr auto VariantBuffer::get() const -> const BuffT & { return this->_buff; }

template <typename T>
bool VariantBuffer::holds_alternative() const noexcept {
  return std::holds_alternative<std::vector<T>>(this->_buff);
}

}  // namespace coolerpp::internal
