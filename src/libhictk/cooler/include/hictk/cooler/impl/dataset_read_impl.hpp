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
  if (offset + num > size()) {
    throw_out_of_range_excp(offset, num);
  }

  buff.resize(num);

  return read(buff.data(), buff.size(), offset);
}

inline std::size_t Dataset::read(std::vector<std::string> &buff, std::size_t num,
                                 std::size_t offset) const {
  if (offset + num > size()) {
    throw_out_of_range_excp(offset, num);
  }

  buff.resize(num);
  const auto str_length = _dataset.getDataType().getSize();
  std::string strbuff(size() * str_length, '\0');
  _dataset.read(strbuff.data(), _dataset.getDataType());

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
    if (offset + num > size()) {
      throw_out_of_range_excp(offset, num);
    }
  }

  using VBuffT = remove_cvref_t<decltype(vbuff.get())>;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using VT = std::variant_alternative_t<i, VBuffT>;
    using T = typename VT::value_type;

    auto read_variant = [&]() {
      if (!vbuff.holds_alternative<T>()) {
        vbuff = std::vector<T>(num);
      }
      return read(vbuff.get<T>(), num, offset);
    };

    auto h5type = get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        return read_variant();
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, num, offset);
    }

    return read_variant();
  }

  unreachable_code();
}
DISABLE_WARNING_POP

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read_n(std::size_t num, std::size_t offset) const {
  BuffT buff{num};
  read(buff, offset);
  return buff;
}

template <typename BuffT, typename T, typename>
inline std::size_t Dataset::read_all(BuffT &buff, std::size_t offset) const {
  const auto num = offset > size() ? std::size_t(0) : size() - offset;
  return read(buff, num, offset);
}

inline hictk::internal::VariantBuffer Dataset::read_all(std::size_t offset) const {
  return read_all<VariantBuffer>(offset);
}

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read_all(std::size_t offset) const {
  BuffT buff{};
  read_all(buff, offset);
  return buff;
}

template <typename N, typename>
inline std::size_t Dataset::read(N &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= size()) {
    throw_out_of_range_excp(offset);
  }

  return read(&buff, 1, offset);
}

inline std::size_t Dataset::read(std::string &buff, std::size_t offset) const {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= size()) {
    throw_out_of_range_excp(offset);
  }

  auto h5type = get_h5type();
  const auto str_length = h5type.getSize();
  buff.resize(str_length);
  select(offset, 1).read(buff.data(), h5type);

  const auto i1 = (std::min)(buff.find('\0'), buff.size());
  buff.resize(i1);

  return offset + 1;
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
template <std::size_t i>
inline std::size_t Dataset::read(GenericVariant &vbuff, std::size_t offset) const {
  if constexpr (i == 0) {
    if (offset >= size()) {
      throw_out_of_range_excp(offset);
    }
  }

  using VBuffT = GenericVariant;
  if constexpr (i < std::variant_size_v<VBuffT>) {
    using T = std::variant_alternative_t<i, VBuffT>;

    auto read_variant = [&]() {
      if (!std::holds_alternative<T>(vbuff)) {
        vbuff = T{};
      }
      return read(std::get<T>(vbuff), offset);
    };

    const auto h5type = get_h5type();
    if constexpr (is_string_v<T>) {
      if (h5type.isFixedLenStr() || h5type.isVariableStr()) {
        return read_variant();
      }
    }
    if (h5type != HighFive::create_datatype<T>()) {
      return read<i + 1>(vbuff, offset);
    }

    return read_variant();
  }

  unreachable_code();
}
DISABLE_WARNING_POP

template <typename BuffT, typename T, typename>
inline BuffT Dataset::read(std::size_t offset) const {
  BuffT buff{};
  read(buff, offset);
  return buff;
}

inline hictk::internal::GenericVariant Dataset::read(std::size_t offset) const {
  return read<GenericVariant>(offset);
}

template <typename BuffT>
inline BuffT Dataset::read_last() const {
  if (empty()) {
    throw_out_of_range_excp(0);
  }
  BuffT buff{};
  read(buff, size() - 1);

  return buff;
}

inline hictk::internal::GenericVariant Dataset::read_last() const {
  return read_last<GenericVariant>();
}

template <typename T>
inline T Dataset::read_attribute(std::string_view key) const {
  return Attribute::read<T>(_dataset, key);
}

template <typename T>
inline void Dataset::read_attribute(std::string_view key, std::vector<T> &buff) const {
  Attribute::read_vector(_dataset, key, buff);
}

inline auto Dataset::read_attribute(std::string_view key, bool missing_ok) const
    -> Attribute::AttributeVar {
  return Attribute::read(_dataset, key, missing_ok);
}

template <typename T>
inline std::size_t Dataset::read(T *buffer, std::size_t buff_size, std::size_t offset) const {
  select(offset, buff_size).read(buffer, HighFive::create_datatype<T>());
  return offset + buff_size;
}

}  // namespace hictk::cooler
