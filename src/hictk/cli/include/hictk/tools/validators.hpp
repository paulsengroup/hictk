// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <exception>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "hictk/cooler/validation.hpp"
#include "hictk/genomic_units.hpp"
#include "hictk/hic/validation.hpp"

namespace hictk::tools {
namespace internal {
class CoolerFileValidator : public CLI::Validator {
 public:
  CoolerFileValidator() : Validator(".[ms]cool") {
    func_ = [](std::string& uri) -> std::string {
      if (cooler::utils::is_cooler(uri) || cooler::utils::is_multires_file(uri) ||
          cooler::utils::is_scool_file(uri)) {
        return "";
      }

      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      return "Not a valid Cooler: " + uri;
    };
  }
};

class SingleResCoolerFileValidator : public CLI::Validator {
 public:
  SingleResCoolerFileValidator() : Validator(".cool") {
    func_ = [](std::string& uri) -> std::string {
      if (!cooler::utils::is_cooler(uri)) {
        if (cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        if (cooler::utils::is_scool_file(uri)) {
          return "URI points to a .scool file: " + uri;
        }
        const auto path = cooler::parse_cooler_uri(uri).file_path;
        if (!std::filesystem::exists(path)) {
          return "No such file: " + path;
        }
        return "Not a valid Cooler: " + uri;
      }
      return "";
    };
  }
};

class MultiresCoolerFileValidator : public CLI::Validator {
 public:
  MultiresCoolerFileValidator() : Validator(".mcool") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!cooler::utils::is_multires_file(uri)) {
        return "Not a valid multi-resolution cooler: " + uri;
      }
      return "";
    };
  }
};

class SingleCellCoolerFileValidator : public CLI::Validator {
 public:
  SingleCellCoolerFileValidator() : Validator(".scool") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!cooler::utils::is_scool_file(uri)) {
        return "Not a valid single-cell cooler: " + uri;
      }
      return "";
    };
  }
};

class HiCFileValidator : public CLI::Validator {
 public:
  HiCFileValidator() : Validator(".hic") {
    func_ = [](std::string& uri) -> std::string {
      const auto path = cooler::parse_cooler_uri(uri).file_path;
      if (!std::filesystem::exists(path)) {
        return "No such file: " + path;
      }
      if (!hic::utils::is_hic_file(path)) {
        return "Not a valid .hic file: " + path;
      }
      return "";
    };
  }
};

class AsGenomicDistanceTransformer : public CLI::Validator {
 public:
  explicit AsGenomicDistanceTransformer() {
    func_ = [](std::string& input) -> std::string {
      try {
        CLI::detail::rtrim(input);
        input = fmt::to_string(parse_genomic_distance<std::uint32_t>(input));
      } catch (const std::exception& e) {
        throw CLI::ValidationError(e.what());
      }
      return {};
    };
  }
};

template <typename Enum>
class StringToEnumChecked : public CLI::Validator {
  static_assert(std::is_enum_v<Enum>);
  [[nodiscard]] static std::string to_lower(std::string_view s) {
    std::string s_lower{s};
    std::transform(s_lower.begin(), s_lower.end(), s_lower.begin(),
                   [](auto c) { return std::tolower(c); });
    return s_lower;
  }

 public:
  // NOLINTNEXTLINE(*-unnecessary-value-param)
  explicit StringToEnumChecked(std::vector<std::pair<std::string, Enum>> mappings) {
    assert(!mappings.empty());

    auto description_formatter = [mappings]() {
      std::vector<std::string> keys(mappings.size());
      std::transform(mappings.begin(), mappings.end(), keys.begin(),
                     [](const auto& kv) { return kv.first; });
      return fmt::format(FMT_STRING("{{{}}}"), fmt::join(keys, ","));
    };

    desc_function_ = description_formatter;

    func_ = [mappings, description_formatter](std::string& input) -> std::string {
      const auto input_lower = to_lower(input);

      const auto match = std::find_if(mappings.begin(), mappings.end(), [&](const auto& kv) {
        return to_lower(kv.first) == input_lower;
      });

      if (match != mappings.end()) {
        return "";
      }

      return fmt::format(FMT_STRING("{} not in {}"), input, description_formatter());
    };
  }
};

}  // namespace internal

// NOLINTBEGIN(cert-err58-cpp)
// clang-format off
inline const auto IsValidHiCFile = internal::HiCFileValidator();
inline const auto IsValidCoolerFile = internal::CoolerFileValidator();
inline const auto IsValidSingleResCoolerFile = internal::SingleResCoolerFileValidator();
inline const auto IsValidMultiresCoolerFile =internal::MultiresCoolerFileValidator();
inline const auto IsValidSingleCellCoolerFile = internal::SingleCellCoolerFileValidator();
inline const auto AsGenomicDistance = internal::AsGenomicDistanceTransformer();

inline const auto ParseHiCMatrixType =
  internal::StringToEnumChecked{
    std::vector<std::pair<std::string, hic::MatrixType>>{
      {"observed", hic::MatrixType::observed},
      {"oe", hic::MatrixType::oe},
      {"expected", hic::MatrixType::expected}
    }
  }.description("");

inline const auto ParseHiCMatrixUnit =
  internal::StringToEnumChecked{
    std::vector<std::pair<std::string, hic::MatrixUnit>>{
      {"BP", hic::MatrixUnit::BP},
      {"FRAG", hic::MatrixUnit::FRAG}
    }
  }.description("");
// clang-format on
// NOLINTEND(cert-err58-cpp)

}  // namespace hictk::tools
