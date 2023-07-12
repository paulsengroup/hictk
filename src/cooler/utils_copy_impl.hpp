// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <H5Lpublic.h>
#include <H5Opublic.h>
#include <H5Ppublic.h>

#include <algorithm>
#include <highfive/H5File.hpp>
#include <string_view>

#include "hictk/cooler.hpp"

namespace hictk::cooler::utils {

inline void copy(std::string_view uri1, std::string_view uri2) {
  try {
    if (!utils::is_cooler(uri1)) {
      throw std::runtime_error("input is not a valid Cooler");
    }
    if (std::filesystem::exists(uri2) && utils::is_cooler(uri2)) {
      throw std::runtime_error("destination already contains a Cooler");
    }

    /* Attempt to open an existing HDF5 file first. Need to open the dst file
         before the src file just in case that the dst and src are the same file
    */
    const auto [path1, grp1] = cooler::parse_cooler_uri(uri1);
    const auto [path2, grp2] = cooler::parse_cooler_uri(uri2);
    const HighFive::File fin(path1, HighFive::File::ReadOnly);

    // TODO check --force overwrite
    HighFive::File fout(path2, HighFive::File::OpenOrCreate);
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

    auto copy_h5_object = [&](const auto &src, const auto &dest, const auto &obj) {
      if (H5Ocopy(src.getId(),    /* Source file or group identifier */
                  obj.data(),     /* Name of the source object to be copied */
                  dest.getId(),   /* Destination file or group identifier  */
                  obj.data(),     /* Name of the destination object  */
                  ocpl_id,        /* Object copy property list */
                  lcpl_id) < 0) { /* Link creation property list */
        throw std::runtime_error(
            fmt::format(FMT_STRING("H5Ocopy failed for {}::/{} failed"), src.getPath(), obj));
      }
    };

    auto copy_h5_attribute = [&](const HighFive::Group &src, HighFive::Group dest,
                                 const auto &attr_name) {
      auto attr = src.getAttribute(attr_name);
      auto attr_out = dest.createAttribute(attr_name, attr.getSpace(), attr.getDataType());

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

    for (const auto &obj : fin.getGroup(grp1).listObjectNames()) {
      copy_h5_object(fin.getGroup(grp1), fout.getGroup(grp2), obj);
    }

    for (const auto &attr : fin.getGroup(grp1).listAttributeNames()) {
      copy_h5_attribute(fin.getGroup(grp1), fout.getGroup(grp2), std::string{attr});
    }

    /* close propertis */
    if (H5Pclose(ocpl_id) < 0) {
      throw std::runtime_error("H5Pclose failed");
    }
    if (H5Pclose(lcpl_id) < 0) {
      throw std::runtime_error("H5Pclose failed");
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to copy Cooler from {} to {}: {}"), uri1, uri2, e.what()));
  }
}

}  // namespace hictk::cooler::utils
