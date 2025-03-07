// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>

namespace hictk::fuzzer {

// NOLINTBEGIN(*-avoid-magic-numbers)
struct Config {
  std::uint16_t task_id{0};
  std::filesystem::path exec{};
  std::filesystem::path test_uri{};
  std::filesystem::path reference_uri{};

  std::uint32_t resolution{};
  double _1d_to_2d_query_ratio{0.33};
  double duration{60.0};
  double query_relative_length_avg{0.03};
  double query_relative_length_std{0.5};
  std::uint64_t min_query_length{500'000};
  std::uint64_t max_query_length{10'000'000};
  std::optional<std::uint64_t> diagonal_band_width{};
  std::string query_format{"df"};
  bool join{false};
  std::string normalization{"NONE"};
  std::optional<std::uint64_t> seed{};
  std::size_t nproc{1};
  bool suppress_python_warnings{true};
  std::int16_t verbosity{3};
};
// NOLINTEND(*-avoid-magic-numbers)

}  // namespace hictk::fuzzer
