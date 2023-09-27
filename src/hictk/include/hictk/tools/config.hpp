// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"

namespace hictk::tools {

struct BalanceConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path tmp_dir{std::filesystem::temp_directory_path()};
  std::filesystem::path juicer_tools_jar{};

  std::string mode{"gw"};
  std::size_t masked_diags{2};
  double mad_max{5.0};
  std::size_t min_nnz{10};
  std::size_t min_count{0};
  double tolerance{1.0e-5};
  std::size_t max_iters{500};
  bool rescale_marginals{true};
  std::string name{"weight"};
  bool in_memory{false};
  bool stdout_{false};
  std::uint8_t zstd_compression_lvl{3};
  std::size_t threads{1};
  std::size_t chunk_size{10'000'000};
  std::size_t juicer_tools_xmx{256'000'000};

  std::uint8_t verbosity{4};
  bool force{false};
};

struct ConvertConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path path_to_output{};
  std::filesystem::path tmp_dir{};
  std::filesystem::path juicer_tools_jar{};
  std::string input_format{};
  std::string output_format{};

  std::vector<std::uint32_t> resolutions{};
  std::string genome{};

  std::vector<balancing::Method> normalization_methods{};
  bool fail_if_normalization_method_is_not_avaliable{false};

  std::uint8_t gzip_compression_lvl{6};
  std::size_t processes{2};

  std::size_t juicer_tools_xmx{32'000'000'000};
  std::uint8_t verbosity{4};
  bool force{false};
};

struct DumpConfig {
  std::string uri{};
  std::string format{};

  std::string range1{"all"};
  std::string range2{"all"};
  std::filesystem::path query_file{};

  std::string table{"pixels"};
  bool join{false};
  bool sorted{true};

  std::string normalization{"NONE"};
  std::string weight_type{"infer"};
  hic::MatrixType matrix_type{hic::MatrixType::observed};
  hic::MatrixUnit matrix_unit{hic::MatrixUnit::BP};
  std::uint32_t resolution{};
  std::uint8_t verbosity{2};
  bool force{false};
};

struct LoadConfig {
  std::string uri{};

  std::filesystem::path path_to_chrom_sizes{};
  std::uint32_t bin_size{};
  std::string format{};
  std::string assembly{"unknown"};
  bool count_as_float{false};
  bool assume_sorted{false};
  bool force{false};
  bool validate_pixels{true};

  std::uint8_t verbosity{4};
  std::size_t batch_size{20'000'000};
};

struct MergeConfig {
  std::vector<std::string> input_uris{};
  std::filesystem::path output_uri{};
  std::filesystem::path tmp_dir{};

  std::size_t chunk_size{5'000'000};

  bool force{false};
  std::uint8_t verbosity{4};
};

struct ValidateConfig {
  std::string uri{};
  bool validate_index{false};
  bool quiet{false};
  std::uint8_t verbosity{4};
};

struct ZoomifyConfig {
  std::string input_uri{};
  std::string output_path{};

  std::vector<std::uint32_t> resolutions{};
  bool copy_base_resolution{true};

  bool force{false};
  std::uint8_t verbosity{4};
};

// clang-format off
using Config = std::variant<std::monostate,
                            BalanceConfig,
                            ConvertConfig,
                            DumpConfig,
                            LoadConfig,
                            MergeConfig,
                            ValidateConfig,
                            ZoomifyConfig>;
// clang-format on

}  // namespace hictk::tools
