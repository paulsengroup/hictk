// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <highfive/H5DataType.hpp>
#include <highfive/H5Utility.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"

namespace hictk::cooler {

template <typename N, typename>
inline std::size_t Dataset::read(std::vector<N> &buff, std::size_t num, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset + num > this->size()) {
    this->throw_out_of_range_excp(offset, num);
  }

  auto h5type = this->get_h5type();
  buff.resize(num);
  this->select(offset, num).read(buff.data(), HighFive::create_datatype<N>());

  return offset + num;
}

inline std::size_t Dataset::read(std::vector<std::string> &buff, std::size_t num,
                                 std::size_t offset) const {
  if (offset + num > this->size()) {
    this->throw_out_of_range_excp(offset, num);
  }

  buff.resize(num);
  const auto str_length = this->_dataset.getDataType().getSize();
  std::string strbuff(this->size() * str_length, '\0');
  this->_dataset.read(strbuff.data(), this->_dataset.getDataType());

  for (std::size_t i = 0; i < buff.size(); ++i) {
    const auto i0 = (offset + i) * str_length;
    const auto i1 = (std::min)(strbuff.find('\0', i0), i0 + str_length);

    buff[i] = strbuff.substr(i0, i1 - i0);
  }

  return offset + buff.size();
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
template <std::size_t i>
inline std::size_t Dataset::read(VariantBuffer &vbuff, std::size_t num, std::size_t offset) const {
  if constexpr (i == 0) {
    if (offset + num > this->size()) {
      this->throw_out_of_range_excp(offset, num);
    }
  }

  using VBuffT = remove_cvref_t<decltype(vbuff.get())>;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using VT = std::variant_alternative_t<i, VBuffT>;
    using T = typename VT::value_type;

    auto h5type = this->get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        goto READ_VARIANT;  // NOLINT
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, num, offset);
    }

#if !defined(__clang__)
    // Workaround for buggy -Wunused-label on GCC and MSVC
    goto READ_VARIANT;  // NOLINT
#endif

  READ_VARIANT:
    if (!std::holds_alternative<VT>(vbuff.get())) {
      vbuff = std::vector<T>(num);
    }

    return this->read(vbuff.get<T>(), num, offset);
  }

  unreachable_code();
}
DISABLE_WARNING_POP

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read_n(std::size_t num, std::size_t offset) const {
  BuffT buff{num};
  this->read(buff, offset);
  return buff;
}

template <typename BuffT, typename T, typename>
inline std::size_t Dataset::read_all(BuffT &buff, std::size_t offset) const {
  const auto num = offset > this->size() ? std::size_t(0) : this->size() - offset;
  return this->read(buff, num, offset);
}

inline hictk::internal::VariantBuffer Dataset::read_all(std::size_t offset) const {
  return this->read_all<VariantBuffer>(offset);
}

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read_all(std::size_t offset) const {
  BuffT buff{};
  this->read_all(buff, offset);
  return buff;
}

template <typename N, typename>
inline std::size_t Dataset::read(N &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= this->size()) {
    this->throw_out_of_range_excp(offset);
  }

  this->select(offset, 1).read(&buff, HighFive::create_datatype<N>());
  return offset + 1;
}

inline std::size_t Dataset::read(std::string &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= this->size()) {
    this->throw_out_of_range_excp(offset);
  }

  auto h5type = this->get_h5type();
  const auto str_length = h5type.getSize();
  buff.resize(str_length);
  this->select(offset, 1).read(buff.data(), h5type);

  const auto i1 = (std::min)(buff.find('\0'), buff.size());
  buff.resize(i1);

  return offset + 1;
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
template <std::size_t i>
inline std::size_t Dataset::read(GenericVariant &vbuff, std::size_t offset) const {
  if constexpr (i == 0) {
    if (offset >= this->size()) {
      this->throw_out_of_range_excp(offset);
    }
  }

  using VBuffT = GenericVariant;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using T = std::variant_alternative_t<i, VBuffT>;

    const auto h5type = this->get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        goto READ_VARIANT;  // NOLINT
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, offset);
    }

#if !defined(__clang__)
    // Workaround for buggy -Wunused-label on GCC and MSVC
    goto READ_VARIANT;  // NOLINT
#endif

  READ_VARIANT:
    if (!std::holds_alternative<T>(vbuff)) {
      vbuff = T{};
    }

    return this->read(std::get<T>(vbuff), offset);
  }

  unreachable_code();
}
DISABLE_WARNING_POP

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read(std::size_t offset) const {
  BuffT buff{};
  this->read(buff, offset);
  return buff;
}

inline hictk::internal::GenericVariant Dataset::read(std::size_t offset) const {
  return this->read<GenericVariant>(offset);
}

template <typename BuffT>
inline BuffT Dataset::read_last() const {
  if (this->empty()) {
    this->throw_out_of_range_excp(0);
  }
  BuffT buff{};
  this->read(buff, this->size() - 1);

  return buff;
}

inline hictk::internal::GenericVariant Dataset::read_last() const {
  return this->read_last<GenericVariant>();
}

template <typename T>
inline T Dataset::read_attribute(std::string_view key) const {
  if constexpr (std::is_same_v<T, bool>) {
    return static_cast<bool>(this->_dataset.getAttribute(std::string{key}).read<int>());
  } else {
    return Attribute::read<T>(this->_dataset, key);
  }
}

template <typename T>
inline void Dataset::read_attribute(std::string_view key, std::vector<T> &buff) const {
  Attribute::read_vector(this->_dataset, key, buff);
}

inline auto Dataset::read_attribute(std::string_view key, bool missing_ok) const
    -> Attribute::AttributeVar {
  return Attribute::read(this->_dataset, key, missing_ok);
}

}  // namespace hictk::cooler
