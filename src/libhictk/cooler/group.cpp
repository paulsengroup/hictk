// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/group.hpp"

#include <fmt/format.h>

#include <highfive/H5Group.hpp>
#include <string>
#include <utility>

namespace hictk::cooler {

RootGroup::RootGroup(HighFive::Group grp) noexcept : _group(std::move(grp)) {}
RootGroup::RootGroup(RootGroup&& other) noexcept : _group(std::move(other._group)) {}
RootGroup& RootGroup::operator=(RootGroup&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  _group = std::move(other._group);

  return *this;
}

std::string RootGroup::file_name() const {
  if (!_group.has_value()) {
    return "";
  }
  return _group->getFile().getName();
}

std::string RootGroup::hdf5_path() const {
  if (!_group.has_value()) {
    return "";
  }
  return _group->getPath();
}

std::string RootGroup::uri() const {
  if (!_group.has_value()) {
    return "";
  }
  return fmt::format(FMT_STRING("{}::{}"), file_name(), hdf5_path());
}

Group::Group(RootGroup root_grp, HighFive::Group grp) noexcept
    : _root_group(std::move(root_grp)), _group(std::move(grp)) {}

}  // namespace hictk::cooler
