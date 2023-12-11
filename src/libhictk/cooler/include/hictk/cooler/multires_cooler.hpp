// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler {

struct MultiResAttributes {
  // Mandatory attributes
  std::string format{MCOOL_MAGIC};
  std::uint8_t format_version{2};

  std::optional<std::string> bin_type{"fixed"};

  MultiResAttributes() = default;
  [[nodiscard]] bool operator==(const MultiResAttributes& other) const noexcept;
  [[nodiscard]] bool operator!=(const MultiResAttributes& other) const noexcept;
};

class MultiResFile {
  std::unique_ptr<RootGroup> _root_grp{};
  std::vector<std::uint32_t> _resolutions{};
  MultiResAttributes _attrs{};
  Reference _chroms{};

  MultiResFile(HighFive::File fp, Reference chroms, std::vector<std::uint32_t> resolutions,
               MultiResAttributes attrs);

 public:
  explicit MultiResFile(const std::filesystem::path& path,
                        unsigned int mode = HighFive::File::ReadOnly);
  [[nodiscard]] static MultiResFile create(const std::filesystem::path& path,
                                           const Reference& chroms, bool force_overwrite = false);
  template <typename ResolutionIt>
  [[nodiscard]] static MultiResFile create(const std::filesystem::path& path, const File& base,
                                           ResolutionIt first_res, ResolutionIt last_res,
                                           bool force_overwrite = false);

  [[nodiscard]] constexpr const std::vector<std::uint32_t>& resolutions() const noexcept;
  [[nodiscard]] constexpr const MultiResAttributes& attributes() const noexcept;
  [[nodiscard]] File open(std::uint32_t resolution) const;
  File copy_resolution(const cooler::File& clr);
  template <typename N = DefaultPixelT>
  File create_resolution(std::uint32_t resolution, Attributes attributes = Attributes::init<N>(0));
  RootGroup init_resolution(std::uint32_t resolution);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] std::string path() const;
  [[nodiscard]] auto chromosomes() const noexcept -> const Reference&;

  [[nodiscard]] static std::uint32_t compute_base_resolution(
      const std::vector<std::uint32_t>& resolutions, std::uint32_t target_res);

  static void coarsen(const File& clr1, File& clr2, std::vector<ThinPixel<std::int32_t>>& buffer);

 private:
  [[nodiscard]] static std::vector<std::uint32_t> read_resolutions(const HighFive::File& f);
  [[nodiscard]] static MultiResAttributes read_attributes(const HighFive::File& f);
};

}  // namespace hictk::cooler

#include "./impl/multires_cooler_impl.hpp"  // NOLINT
