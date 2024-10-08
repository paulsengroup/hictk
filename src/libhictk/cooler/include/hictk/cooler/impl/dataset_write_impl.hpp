// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#if __has_include(<hdf5/hdf5.h>)
#include <hdf5/H5Dpublic.h>
#include <hdf5/H5Ipublic.h>
#include <hdf5/H5Ppublic.h>
#include <hdf5/H5Tpublic.h>
#else
#include <H5Dpublic.h>
#include <H5Ipublic.h>
#include <H5Ppublic.h>
#include <H5Tpublic.h>
#endif

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Utility.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictk/cooler/attribute.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

inline std::size_t Dataset::write(const std::vector<std::string> &buff, std::size_t offset,
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

template <typename N, typename>
inline std::size_t Dataset::write(const std::vector<N> &buff, std::size_t offset,
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

  select(offset, buff.size()).write(buff);
  return size();
}

// NOLINTNEXTLINE(*-convert-member-functions-to-static)
inline std::size_t Dataset::write(const VariantBuffer &vbuff, std::size_t offset,
                                  bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit([&](const auto &buff) { new_offset = write(buff, offset, allow_dataset_resize); },
             vbuff.get());

  return new_offset;
}

template <typename InputIt, typename UnaryOperation,  // NOLINTNEXTLINE(*modernize-type-traits)
          typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int> *>
inline std::size_t Dataset::write(InputIt first_value, const InputIt &last_value,
                                  std::size_t offset, bool allow_dataset_resize,
                                  UnaryOperation op) {
  if (first_value == last_value) {
    return offset;
  }
  using T = remove_cvref_t<decltype(op(*first_value))>;
  constexpr std::size_t buffer_capacity =
      is_string_v<T> ? 256 : (1ULL << 20U) / sizeof(std::uint64_t);
  if (_buff.holds_alternative<T>()) {
    _buff.resize<T>(buffer_capacity);
  } else {
    _buff = std::vector<T>(buffer_capacity);
  }

  auto &buff = _buff.get<T>();
  buff.clear();

  while (first_value != last_value) {
    if (buff.size() == buff.capacity()) {
      write(buff, offset, allow_dataset_resize);
      offset += buff.size();
      buff.clear();
    }

    buff.emplace_back(op(*first_value));
    std::ignore = ++first_value;
  }

  if (!buff.empty()) {
    write(buff, offset, allow_dataset_resize);
    offset += buff.size();
  }

  return offset;
}

template <typename InputIt, typename UnaryOperation,  // NOLINTNEXTLINE(*modernize-type-traits)
          typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int> *>
inline std::size_t Dataset::append(InputIt first_value, const InputIt &last_value,
                                   UnaryOperation op) {
  return write(std::move(first_value), last_value, size(), true, op);
}

template <typename N, typename>
inline std::size_t Dataset::write(N buff, std::size_t offset, bool allow_dataset_resize) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (offset >= size()) {
    if (allow_dataset_resize) {
      resize(offset + 1);
    } else {
      throw_out_of_range_excp(offset);
    }
  }

  select(offset).write(buff);
  return ++_dataset_size;
}

inline std::size_t Dataset::write(std::string buff, std::size_t offset, bool allow_dataset_resize) {
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
inline std::size_t Dataset::write(const GenericVariant &vbuff, std::size_t offset,
                                  bool allow_dataset_resize) {
  std::size_t new_offset{};
  std::visit([&](const auto &buff) { new_offset = write(buff, offset, allow_dataset_resize); },
             vbuff);

  return new_offset;
}

template <typename BuffT>
inline std::size_t Dataset::append(const BuffT &buff) {
  return write(buff, size(), true);
}

template <typename T>
inline void Dataset::write_attribute(std::string_view key, const T &value,
                                     bool overwrite_if_exists) {
  Attribute::write(_dataset, key, value, overwrite_if_exists);
}

inline HighFive::DataSet Dataset::create_fixed_str_dataset(
    RootGroup &root_grp, std::string_view path, std::size_t max_str_length, std::size_t max_dim,
    const HighFive::DataSetAccessProps &aprops, const HighFive::DataSetCreateProps &cprops) {
  assert(max_str_length != 0);

  const auto [group_name, dataset_name] = parse_uri(path);
  auto group = root_grp().getGroup(group_name);
  if (group.exist(dataset_name)) {
    throw std::runtime_error(fmt::format(FMT_STRING("Dataset at URI \"{}\" already exists"), path));
  }

  auto dspace = HighFive::DataSpace({0}, {max_dim});

  // Unfortunately we have to drop down to the C api to create this kind of dataset for the time
  // being
  auto dtype_id = H5Tcopy(H5T_C_S1);
  H5Tset_cset(dtype_id, H5T_CSET_ASCII);
  H5Tset_size(dtype_id, max_str_length);
  H5Tset_strpad(dtype_id, H5T_STR_NULLPAD);

  const auto hid = H5Dcreate(group.getId(), dataset_name.c_str(), dtype_id, dspace.getId(),
                             H5P_DEFAULT, cprops.getId(), aprops.getId());

  if (hid == H5I_INVALID_HID) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to create dataset at URI \"{}\""), path));
  }

  std::ignore = H5Dclose(hid);
  return group.getDataSet(dataset_name);
}

}  // namespace hictk::cooler
