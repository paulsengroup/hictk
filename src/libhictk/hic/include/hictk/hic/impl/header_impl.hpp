// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstddef>
#include <functional>
#include <ios>
#include <stdexcept>

#include "hictk/hash.hpp"
#include "hictk/hic/filestream.hpp"

namespace hictk::hic::internal {

constexpr HiCHeader::operator bool() const noexcept { return footerPosition >= 0; }

inline bool HiCHeader::operator==(const HiCHeader &other) const noexcept {
  return url == other.url && footerPosition == other.footerPosition;
}

inline bool HiCHeader::operator!=(const HiCHeader &other) const noexcept {
  return !(*this == other);
}

inline std::string HiCHeader::serialize(BinaryBuffer &buffer, bool clear) const {
  if (version != 9) {
    throw std::runtime_error("serializing header for file version other than v9 is not supported.");
  }
  if (chromosomes.empty()) {
    throw std::runtime_error("serializing a header without chromosomes is not supported.");
  }

  if (clear) {
    buffer.clear();
  }

  buffer.write("HIC\0", 4);
  buffer.write(version);
  buffer.write(footerPosition);
  buffer.write(genomeID.c_str(), genomeID.size() + 1);
  buffer.write(normVectorIndexPosition);
  buffer.write(normVectorIndexLength);

  // Write attributes
  const auto nAttributes = static_cast<std::int32_t>(attributes.size());
  buffer.write(nAttributes);
  for (const auto &[k, v] : attributes) {
    buffer.write(k.c_str(), k.size() + 1);
    buffer.write(v.c_str(), v.size() + 1);
  }

  // Write chromosomes
  auto numChromosomes = static_cast<std::int32_t>(chromosomes.size());
  buffer.write(numChromosomes);

  for (const Chromosome &c : chromosomes) {
    const auto name = std::string{c.name()};
    buffer.write(name.c_str(), name.size() + 1);
    buffer.write<std::int64_t>(c.size());
  }

  // write resolutions
  buffer.write(static_cast<std::int32_t>(resolutions.size()));
  const std::vector<std::int32_t> resolutions_(resolutions.begin(), resolutions.end());
  buffer.write(resolutions_);

  // write fragments: TODO
  const std::int32_t nFragResolutions = 0;
  buffer.write(nFragResolutions);

  return buffer.get();
}

inline HiCHeader HiCHeader::deserialize(filestream::FileStream &fs) {
  fs.seekg(0, std::ios::beg);
  const auto magic_string_found = fs.getline('\0') == "HIC";
  if (!magic_string_found) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Hi-C magic string is missing. {} does not appear to be a hic file"),
                    fs.path()));
  }

  HiCHeader header{fs.path()};

  fs.read(header.version);
  if (header.version < 6) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(".hic version 5 and older are no longer supported. Found version {}"),
        header.version));
  }
  fs.read(header.footerPosition);
  if (header.footerPosition < 0 || header.footerPosition >= static_cast<std::int64_t>(fs.size())) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file appears to be corrupted: expected footerPosition to be "
                               "between 0 and {}, found {}"),
                    fs.size(), header.footerPosition));
  }

  fs.getline(header.genomeID, '\0');
  if (header.genomeID.empty()) {
    header.genomeID = "unknown";
  }

  if (header.version > 8) {
    fs.read(header.normVectorIndexPosition);
    fs.read(header.normVectorIndexLength);
  }

  const auto nAttributes = fs.read<std::int32_t>();

  // reading attribute-value dictionary
  for (std::int32_t i = 0; i < nAttributes; i++) {
    auto key = fs.getline('\0');    // key
    auto value = fs.getline('\0');  // value
    header.attributes.emplace(std::move(key), std::move(value));
  }

  // Read chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(fs.read<std::int32_t>());
  std::vector<std::string> chrom_names(numChromosomes);
  std::vector<std::uint32_t> chrom_sizes(numChromosomes);
  for (std::size_t i = 0; i < chrom_names.size(); ++i) {
    fs.getline(chrom_names[i], '\0');
    chrom_sizes[i] = static_cast<std::uint32_t>(
        header.version > 8 ? fs.read<std::int64_t>()
                           : static_cast<std::int64_t>(fs.read<std::int32_t>()));
  }

  if (chrom_names.empty()) {
    throw std::runtime_error("unable to read chromosomes");
  }

  header.chromosomes = Reference(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());

  // Read resolutions
  const auto numResolutions = static_cast<std::size_t>(fs.read<std::int32_t>());
  if (numResolutions == 0) {
    throw std::runtime_error("unable to read the list of available resolutions");
  }
  header.resolutions.resize(numResolutions);
  std::generate(header.resolutions.begin(), header.resolutions.end(), [&]() {
    const auto res = fs.read<std::int32_t>();
    assert(res > 0);
    return static_cast<std::uint32_t>(res);
  });

  return header;
}

}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::HiCHeader> {
  inline std::size_t operator()(hictk::hic::internal::HiCHeader const &h) const noexcept {
    return hictk::internal::hash_combine(0, h.url, h.footerPosition);
  }
};
