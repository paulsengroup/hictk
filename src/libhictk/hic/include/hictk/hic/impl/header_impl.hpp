// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ios>
#include <stdexcept>

#include "hictk/filestream.hpp"
#include "hictk/hash.hpp"

namespace hictk::hic::internal {

constexpr HiCHeader::operator bool() const noexcept { return footerPosition >= 0; }

inline bool HiCHeader::operator==(const HiCHeader &other) const noexcept {
  return url == other.url && footerPosition == other.footerPosition;
}

inline bool HiCHeader::operator!=(const HiCHeader &other) const noexcept {
  return !(*this == other);
}

inline std::string HiCHeader::serialize(BinaryBuffer &buffer, bool clear) const {
  if (version != 9) {  // NOLINT(*-avoid-magic-numbers)
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
  buffer.write(genomeID, true);
  buffer.write(normVectorIndexPosition);
  buffer.write(normVectorIndexLength);

  // Write attributes
  const auto nAttributes = static_cast<std::int32_t>(attributes.size());
  buffer.write(nAttributes);
  for (const auto &[k, v] : attributes) {
    buffer.write(k, true);
    buffer.write(v, true);
  }

  // Write chromosomes
  auto numChromosomes = static_cast<std::int32_t>(chromosomes.size());
  buffer.write(numChromosomes);

  for (const Chromosome &c : chromosomes) {
    const auto name = std::string{c.name()};
    buffer.write(name, true);
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

inline HiCHeader HiCHeader::deserialize(std::streampos offset, filestream::FileStream<> &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

inline HiCHeader HiCHeader::unsafe_deserialize(std::streampos offset,
                                               filestream::FileStream<> &fs) {
  fs.unsafe_seekg(offset);
  std::string strbuff;
  fs.unsafe_getline(strbuff, '\0');
  const auto magic_string_found = strbuff == "HIC";
  if (!magic_string_found) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Hi-C magic string is missing. {} does not appear to be a hic file"),
                    fs.path()));
  }

  HiCHeader header{fs.path()};

  fs.unsafe_read(header.version);
  if (header.version < 6) {  // NOLINT(*-avoid-magic-numbers)
    throw std::runtime_error(fmt::format(FMT_STRING("unable to open .hic file with version={}: "
                                                    "version 5 and older are no longer supported"),
                                         header.version));
  }
  if (header.version > 9) {  // NOLINT(*-avoid-magic-numbers)
    throw std::runtime_error(fmt::format(FMT_STRING("unable to open .hic file with version={}: "
                                                    "versions newer than v9 are not yet supported"),
                                         header.version));
  }
  fs.unsafe_read(header.footerPosition);
  if (header.footerPosition < 0 ||
      header.footerPosition >= conditional_static_cast<std::int64_t>(fs.unsafe_size())) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file appears to be corrupted: expected footerPosition to be "
                               "between 0 and {}, found {}"),
                    fs.unsafe_size(), header.footerPosition));
  }

  fs.unsafe_getline(header.genomeID, '\0');
  if (header.genomeID.empty()) {
    header.genomeID = "unknown";
  }

  if (header.version > 8) {  // NOLINT(*-avoid-magic-numbers)
    fs.unsafe_read(header.normVectorIndexPosition);
    fs.unsafe_read(header.normVectorIndexLength);
  }

  const auto nAttributes = fs.unsafe_read<std::int32_t>();

  // reading attribute-value dictionary
  std::string key;
  std::string value;
  for (std::int32_t i = 0; i < nAttributes; i++) {
    fs.unsafe_getline(key, '\0');
    fs.unsafe_getline(value, '\0');
    header.attributes.emplace(key, value);
  }

  // Read chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(fs.unsafe_read<std::int32_t>());
  std::vector<std::string> chrom_names(numChromosomes);
  std::vector<std::uint32_t> chrom_sizes(numChromosomes);
  for (std::size_t i = 0; i < chrom_names.size(); ++i) {
    fs.unsafe_getline(chrom_names[i], '\0');
    chrom_sizes[i] = static_cast<std::uint32_t>(
        header.version > 8 ? fs.unsafe_read<std::int64_t>()  // NOLINT(*-avoid-magic-numbers)
                           : static_cast<std::int64_t>(fs.unsafe_read<std::int32_t>()));
  }

  if (chrom_names.empty()) {
    throw std::runtime_error("unable to read chromosomes");
  }

  header.chromosomes = Reference(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());

  // Read resolutions
  const auto numResolutions = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
  if (numResolutions == 0) {
    throw std::runtime_error("unable to read the list of available resolutions");
  }

  // sometimes .hic files have duplicate resolutions for some obscure reason...
  phmap::btree_set<std::uint32_t> resolutions{};
  for (std::size_t i = 0; i < numResolutions; ++i) {
    const auto res = static_cast<std::uint32_t>(fs.unsafe_read<std::int32_t>());
    resolutions.emplace(res);
  }

  header.resolutions.reserve(resolutions.size());
  std::copy(resolutions.begin(), resolutions.end(), std::back_inserter(header.resolutions));

  return header;
}

}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::HiCHeader> {
  std::size_t operator()(hictk::hic::internal::HiCHeader const &h) const noexcept {
    return hictk::internal::hash_combine(0, h.url, h.footerPosition);
  }
};
