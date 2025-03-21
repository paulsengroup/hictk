// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/hic/common.hpp"

namespace hictk::tools {

// NOLINTBEGIN(*-avoid-magic-numbers)

static constexpr std::int16_t DEFAULT_COOL_COMPRESSION_LEVEL = 6;
static constexpr std::int16_t MAX_COOL_COMPRESSION_LEVEL = 9;

static constexpr std::int16_t DEFAULT_HIC_COMPRESSION_LEVEL = 10;
static constexpr std::int16_t MAX_HIC_COMPRESSION_LEVEL = 12;

static constexpr std::int16_t DEFAULT_ZSTD_COMPRESSION_LEVEL = 3;
static constexpr std::int16_t MAX_ZSTD_COMPRESSION_LEVEL = 19;

struct BalanceICEConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path tmp_dir{};

  std::string mode{"gw"};
  std::size_t masked_diags{2};
  double mad_max{5.0};
  std::size_t min_nnz{10};
  std::size_t min_count{0};
  double tolerance{1.0e-5};
  std::size_t max_iters{500};
  bool rescale_marginals{true};
  std::string name{};
  bool in_memory{false};
  bool symlink_to_weight{true};
  bool stdout_{false};
  std::int16_t zstd_compression_lvl{DEFAULT_ZSTD_COMPRESSION_LEVEL};
  std::size_t threads{1};
  std::size_t chunk_size{10'000'000};

  std::int16_t verbosity{3};
  bool force{false};
};

struct BalanceSCALEConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path tmp_dir{};

  std::string mode{"gw"};
  double max_percentile{10};
  double max_row_sum_error{0.05};
  double tolerance{1.0e-4};
  std::size_t max_iters{500};
  bool rescale_marginals{true};
  std::string name{};
  bool in_memory{false};
  bool symlink_to_weight{true};
  bool stdout_{false};
  std::int16_t zstd_compression_lvl{DEFAULT_ZSTD_COMPRESSION_LEVEL};
  std::size_t threads{1};
  std::size_t chunk_size{10'000'000};

  std::int16_t verbosity{3};
  bool force{false};
};

struct BalanceVCConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path tmp_dir{};  // unused

  std::string mode{"gw"};
  bool rescale_marginals{true};
  std::string name{};
  bool symlink_to_weight{true};
  bool stdout_{false};

  std::int16_t verbosity{3};
  bool force{false};
};

struct ConvertConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path path_to_output{};
  std::filesystem::path tmp_dir{};
  std::string input_format{};
  std::string output_format{};
  std::string count_type{"auto"};

  std::vector<std::uint32_t> resolutions{};
  std::string genome{};

  std::vector<balancing::Method> normalization_methods{};
  bool fail_if_normalization_method_is_not_avaliable{false};
  bool skip_all_vs_all_matrix{false};

  std::uint32_t compression_lvl{DEFAULT_COOL_COMPRESSION_LEVEL};
  std::size_t threads{2};
  std::size_t chunk_size{10'000'000};

  std::int16_t verbosity{3};
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

  bool cis_only{false};
  bool trans_only{false};

  std::string normalization{"NONE"};
  hic::MatrixType matrix_type{hic::MatrixType::observed};
  hic::MatrixUnit matrix_unit{hic::MatrixUnit::BP};
  std::optional<std::uint32_t> resolution{};
  std::int16_t verbosity{2};
  bool force{false};
};

struct FixMcoolConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path path_to_output{};
  std::filesystem::path tmp_dir{};

  bool skip_balancing{false};
  bool check_base_resolution{false};

  bool in_memory{false};
  std::int16_t zstd_compression_lvl{DEFAULT_ZSTD_COMPRESSION_LEVEL};
  std::size_t chunk_size{10'000'000};

  std::size_t threads{1};
  std::int16_t verbosity{3};
  bool force{false};
};

struct LoadConfig {
  std::filesystem::path input_path{"-"};
  std::string output_path{};

  std::filesystem::path path_to_chrom_sizes{};
  std::filesystem::path path_to_bin_table{};
  std::filesystem::path tmp_dir{};
  std::uint32_t bin_size{};

  std::string format{};
  std::string assembly{"unknown"};
  bool drop_unknown_chroms{false};
  bool one_based{true};
  std::int64_t offset{0};
  bool count_as_float{false};
  bool assume_sorted{false};
  bool force{false};
  bool validate_pixels{true};
  bool transpose_lower_triangular_pixels{false};
  bool skip_all_vs_all_matrix{true};

  std::string output_format{"auto"};

  std::size_t threads{2};
  std::uint32_t compression_lvl{DEFAULT_COOL_COMPRESSION_LEVEL};

  std::int16_t verbosity{3};
  std::size_t batch_size{10'000'000};
};

struct MergeConfig {
  std::vector<std::string> input_files{};
  std::filesystem::path output_file{};
  std::string output_format{};
  std::optional<std::uint32_t> resolution{};

  std::filesystem::path tmp_dir{};

  std::size_t chunk_size{10'000'000};
  std::uint32_t compression_lvl{DEFAULT_COOL_COMPRESSION_LEVEL};
  std::size_t threads{1};
  bool skip_all_vs_all_matrix{true};
  std::string count_type{"int"};

  bool force{false};
  std::int16_t verbosity{3};
};

struct MetadataConfig {
  std::filesystem::path uri{};
  std::string input_format{};
  std::string output_format{"json"};
  bool include_file_path{false};
  bool recursive{false};

  std::int16_t verbosity{2};
};

struct RenameChromosomesConfig {
  std::string uri{};
  std::filesystem::path path_to_name_mappings{};
  bool add_chr_prefix{false};
  bool remove_chr_prefix{false};
  std::int16_t verbosity{3};
};

struct ValidateConfig {
  std::string uri{};
  bool validate_index{false};
  bool validate_pixels{false};
  std::string output_format{"json"};
  bool include_file_path{true};
  bool exhaustive{true};
  bool quiet{false};
  std::int16_t verbosity{3};
};

struct ZoomifyConfig {
  std::filesystem::path path_to_input{};
  std::filesystem::path path_to_output{};
  std::string input_format{};
  std::string output_format{};
  std::filesystem::path tmp_dir{};

  std::vector<std::uint32_t> resolutions{};
  bool copy_base_resolution{true};
  bool nice_resolution_steps{true};

  std::uint32_t compression_lvl{DEFAULT_COOL_COMPRESSION_LEVEL};
  std::uint32_t threads{1};
  std::size_t batch_size{10'000'000};
  bool skip_all_vs_all_matrix{false};

  bool force{false};
  std::int16_t verbosity{3};
};

// NOLINTEND(*-avoid-magic-numbers)

// clang-format off
using Config = std::variant<std::monostate,
                            BalanceICEConfig,
                            BalanceSCALEConfig,
                            BalanceVCConfig,
                            ConvertConfig,
                            DumpConfig,
                            FixMcoolConfig,
                            LoadConfig,
                            MergeConfig,
                            MetadataConfig,
                            RenameChromosomesConfig,
                            ValidateConfig,
                            ZoomifyConfig>;
// clang-format on

}  // namespace hictk::tools
