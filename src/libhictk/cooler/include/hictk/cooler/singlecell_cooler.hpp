// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <memory>
#include <optional>
#include <string>
#include <string_view>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler {

struct SingleCellAttributes {
  friend File;

  // Mandatory attributes
  std::uint32_t bin_size{0};
  std::string bin_type{"fixed"};
  std::string format{SCOOL_MAGIC};
  std::uint8_t format_version{1};

  // Reserved attributes
  std::optional<std::string> creation_date{Attributes::generate_creation_date()};
  std::optional<std::string> generated_by{HICTK_VERSION_STRING_LONG};
  std::optional<std::string> assembly{"unknown"};
  std::optional<std::string> metadata{"{}"};

  // Optional but common
  std::optional<std::string> format_url{"https://github.com/open2c/cooler"};
  std::optional<std::int64_t> nbins{0};
  std::optional<std::int32_t> ncells{0};
  std::optional<std::int32_t> nchroms{0};
  std::optional<std::string> storage_mode{"symmetric-upper"};

  [[nodiscard]] static SingleCellAttributes init(std::uint32_t bin_size_);
  [[nodiscard]] static SingleCellAttributes init_empty() noexcept;
  [[nodiscard]] bool operator==(const SingleCellAttributes& other) const noexcept;
  [[nodiscard]] bool operator!=(const SingleCellAttributes& other) const noexcept;

 private:
  // Use the init factory methods to construct an Attribute object
  SingleCellAttributes() = default;
};

class SingleCellFile {
  std::unique_ptr<RootGroup> _root_grp{};
  phmap::btree_set<std::string> _cells{};
  SingleCellAttributes _attrs{};
  std::shared_ptr<const BinTable> _bins{};

  SingleCellFile(HighFive::File fp, BinTable bins, SingleCellAttributes attrs);

 public:
  explicit SingleCellFile(const std::filesystem::path& path,
                          unsigned int mode = HighFive::File::ReadOnly);
  [[nodiscard]] static SingleCellFile create(const std::filesystem::path& path,
                                             const Reference& chroms, std::uint32_t bin_size,
                                             bool force_overwrite);
  [[nodiscard]] static SingleCellFile create(const std::filesystem::path& path, BinTable bins,
                                             bool force_overwrite = false);

  [[nodiscard]] constexpr const phmap::btree_set<std::string>& cells() const noexcept;
  [[nodiscard]] constexpr const SingleCellAttributes& attributes() const noexcept;
  [[nodiscard]] File open(std::string_view cell) const;
  template <typename N>
  File create_cell(std::string_view cell, Attributes attrs = Attributes::init<N>(0));

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] std::string path() const;
  [[nodiscard]] auto chromosomes() const noexcept -> const Reference&;
  [[nodiscard]] auto bins() const noexcept -> const BinTable&;
  [[nodiscard]] std::uint32_t bin_size() const noexcept;

  template <typename N>  // NOLINTNEXTLINE(*-use-nodiscard)
  File aggregate(std::string_view uri, bool overwrite_if_exists = false,
                 std::size_t chunk_size = 500'000, std::size_t update_frequency = 10'000'000) const;

 private:
  [[nodiscard]] static SingleCellAttributes read_standard_attributes(const HighFive::File& f,
                                                                     bool initialize_missing);
  [[nodiscard]] static BinTable init_bin_table(const HighFive::File& f);
  [[nodiscard]] static phmap::btree_set<std::string> read_cells(const HighFive::File& f);

  static void create_groups(RootGroup& root_grp);
  static void create_datasets(RootGroup& root_grp, const BinTable& bins);
  static void write_standard_attributes(RootGroup& root_grp, const SingleCellAttributes& attrs);

  template <typename PixelT>
  static void create_cell_datasets(RootGroup& root_grp, std::size_t cache_size_bytes, double w0);
};

}  // namespace hictk::cooler

#include "./impl/singlecell_cooler_impl.hpp"  // NOLINT
