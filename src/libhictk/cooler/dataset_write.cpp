// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cstddef>
#include <highfive/H5Utility.hpp>
#include <string>
#include <variant>
#include <vector>

#include "hictk/cooler/dataset.hpp"

namespace hictk::cooler {

std::size_t Dataset::write(const std::vector<std::string> &buff, std::size_t offset,
                           bool allow_dataset_resize) {
  if (buff.empty()) {
    return offset;
  }
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset + buff.size() > size()) {
    if (allow_dataset_resize) {
      resize(offset + buff.size());
    } else {
      throw_out_of_range_excp(offset, buff.size());
    }
  }

  const auto str_length = get_h5type().getSize();
  std::string strbuff(str_length * buff.size(), '\0');
  for (std::size_t i = 0; i < buff.size(); ++i) {
    strbuff.insert(i * str_length, buff[i]);
  }
  auto dspace = select(offset, buff.size());
  dspace.write_raw(strbuff.data(), dspace.getDataType());

  return size();
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
std::size_t Dataset::write(const VariantBuffer &vbuff, std::size_t offset,
                           bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit([&](const auto &buff) { new_offset = write(buff, offset, allow_dataset_resize); },
             vbuff.get());

  return new_offset;
}

std::size_t Dataset::write(std::string buff, std::size_t offset, bool allow_dataset_resize) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= size()) {
    if (allow_dataset_resize) {
      resize(offset + 1);
    } else {
      throw_out_of_range_excp(offset);
    }
  }

  auto selector = select(offset);

  const auto str_length = selector.getDataType().getSize();
  buff.resize(str_length, '\0');

  selector.write_raw(buff.data(), selector.getDataType());

  return ++_dataset_size;
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
std::size_t Dataset::write(const GenericVariant &vbuff, std::size_t offset,
                           bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit([&](const auto &buff) { new_offset = write(buff, offset, allow_dataset_resize); },
             vbuff);

  return new_offset;
}

}  // namespace hictk::cooler
