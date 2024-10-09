// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#if __has_include(<hdf5/hdf5.h>)
#include <hdf5/H5Opublic.h>
#include <hdf5/H5Ppublic.h>
#else
#include <H5Opublic.h>
#include <H5Ppublic.h>
#endif

#include <exception>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <stdexcept>
#include <string>
#include <string_view>

#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"

namespace hictk::cooler::utils {

inline void copy(std::string_view uri1, std::string_view uri2) {
  const auto uri = cooler::parse_cooler_uri(uri2);
  if (std::filesystem::exists(uri2) && utils::is_cooler(uri2)) {
    throw std::runtime_error("destination already contains a Cooler");
  }
  const HighFive::File dest(uri.file_path, HighFive::File::OpenOrCreate);
  return copy(uri1, RootGroup{dest.getGroup(uri.group_path)});
}

inline void copy(std::string_view uri1, RootGroup dest) {
  try {
    if (!utils::is_cooler(uri1)) {
      throw std::runtime_error("input is not a valid Cooler");
    }

    /* Attempt to open an existing HDF5 file first. Need to open the dst file
         before the src file just in case that the dst and src are the same file
    */
    const auto uri = cooler::parse_cooler_uri(uri1);
    const HighFive::File fin(uri.file_path, HighFive::File::ReadOnly);

    /* create property to pass copy options */
    auto ocpl_id = H5Pcreate(H5P_OBJECT_COPY);
    if (ocpl_id < 0) {
      throw std::runtime_error("H5Pcreate failed");
    }

    /* Create link creation property list */
    auto lcpl_id = H5Pcreate(H5P_LINK_CREATE);
    if (lcpl_id < 0) {
      throw std::runtime_error("Could not create link creation property list: H5Pcreate failed");
    }

    /* Set the intermediate group creation property */
    if (H5Pset_create_intermediate_group(lcpl_id, 1) < 0) {
      throw std::runtime_error(
          "Could not set property for creating parent groups: H5Pset_create_intermediate_group "
          "failed");
    }

    auto copy_h5_object = [&](const auto &src, const auto &dest_, const auto &obj) {
      if (H5Ocopy(src.getId(),    /* Source file or group identifier */
                  obj.data(),     /* Name of the source object to be copied */
                  dest_.getId(),  /* Destination file or group identifier  */
                  obj.data(),     /* Name of the destination object  */
                  ocpl_id,        /* Object copy property list */
                  lcpl_id) < 0) { /* Link creation property list */
        throw std::runtime_error(
            fmt::format(FMT_STRING("H5Ocopy failed for {}::/{} failed"), src.getPath(), obj));
      }
    };

    auto copy_h5_attribute = [&](const HighFive::Group &src, HighFive::Group dest_,
                                 const auto &attr_name) {
      auto attr = src.getAttribute(attr_name);
      auto attr_out = dest_.createAttribute(attr_name, attr.getSpace(), attr.getDataType());

      std::string buffer(attr.getStorageSize(), '\0');
      if (attr.getDataType().isFixedLenStr() || attr.getDataType().isVariableStr()) {
        // Let HighFive deal with strings
        attr.read(buffer);
        attr_out.write(buffer);
      } else {
        attr.read(buffer.data(), attr.getDataType());
        attr_out.write_raw(buffer.data(), attr.getDataType());
      }
    };

    for (const auto &obj : fin.getGroup(uri.group_path).listObjectNames()) {
      copy_h5_object(fin.getGroup(uri.group_path), dest(), obj);
    }

    for (const auto &attr : fin.getGroup(uri.group_path).listAttributeNames()) {
      copy_h5_attribute(fin.getGroup(uri.group_path), dest(), std::string{attr});
    }

    /* close propertis */
    if (H5Pclose(ocpl_id) < 0) {
      throw std::runtime_error("H5Pclose failed");
    }
    if (H5Pclose(lcpl_id) < 0) {
      throw std::runtime_error("H5Pclose failed");
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to copy Cooler from {} to {}: {}"),
                                         uri1, dest.uri(), e.what()));
  }
}

}  // namespace hictk::cooler::utils
