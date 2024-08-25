// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>

namespace hictk::fuzzer {

struct Config {
  std::filesystem::path test_uri{};
  std::filesystem::path reference_uri{};

  std::uint32_t resolution{};
  double _1d_to_2d_query_ratio{0.33};
  double duration{60.0};
  double query_length_avg{1.0e6};
  double query_length_std{250.0e3};
  std::string query_format{"df"};
  bool join{false};
  std::string normalization{"NONE"};
  std::uint64_t seed{11261741397133096960ULL};
  std::size_t nproc{1};
};

}  // namespace hictk::fuzzer
