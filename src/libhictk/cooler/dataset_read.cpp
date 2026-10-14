// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstddef>
#include <highfive/H5Utility.hpp>
#include <highfive/H5Version.hpp>
#include <string>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/generic_variant.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk::cooler {

std::size_t Dataset::read(std::vector<std::string> &buff, std::size_t num,
                          std::size_t offset) const {
  if (offset + num > size()) {
    throw_out_of_range_excp(offset, num);
  }

  buff.resize(num);
  const auto str_length = _dataset.getDataType().getSize();
  std::string strbuff(size() * str_length, '\0');

#if HIGHFIVE_VERSION_MAJOR > 2
  _dataset.read_raw(strbuff.data(), _dataset.getDataType());
#else
  _dataset.read(strbuff.data(), _dataset.getDataType());
#endif

  for (std::size_t i = 0; i < buff.size(); ++i) {
    const auto i0 = (offset + i) * str_length;
    const auto i1 = (std::min)(strbuff.find('\0', i0), i0 + str_length);

    buff[i] = strbuff.substr(i0, i1 - i0);
  }

  return offset + buff.size();
}

hictk::internal::VariantBuffer Dataset::read_all(std::size_t offset) const {
  return read_all<VariantBuffer>(offset);
}

std::size_t Dataset::read(std::string &buff, std::size_t offset) const {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= size()) {
    throw_out_of_range_excp(offset);
  }

  auto h5type = get_h5type();
  const auto str_length = h5type.getSize();
  buff.resize(str_length);
#if HIGHFIVE_VERSION_MAJOR > 2
  select(offset, 1).read_raw(buff.data(), h5type);
#else
  select(offset, 1).read(buff.data(), h5type);
#endif

  const auto i1 = (std::min)(buff.find('\0'), buff.size());
  buff.resize(i1);

  return offset + 1;
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
hictk::internal::GenericVariant Dataset::read(std::size_t offset) const {
  return read<GenericVariant>(offset);
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
hictk::internal::GenericVariant Dataset::read_last() const { return read_last<GenericVariant>(); }

auto Dataset::read_attribute(std::string_view key, bool missing_ok) const
    -> Attribute::AttributeVar {
  return Attribute::read(_dataset, key, missing_ok);
}

}  // namespace hictk::cooler
