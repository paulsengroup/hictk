// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./validate.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstdio>
#include <filesystem>
#include <string>
#include <string_view>

#include "hictk/cooler/validation.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/file_attributes_formatting.hpp"
#include "hictk/tools/toml.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

static void print_report(const toml::table& status, std::string_view format) {
  if (format == "json") {
    fmt::print(FMT_STRING("{}\n"), io::toml::format_to_json(status, {}));
    return;
  }

  if (format == "toml") {
    fmt::print(FMT_STRING("{}\n"), io::toml::format_to_toml(status, {}));
    return;
  }

  assert(format == "yaml");
  fmt::print(FMT_STRING("{}\n"), io::toml::format_to_yaml(status, {}));
}

[[nodiscard]] static toml::table merge_tables(const toml::table& t1, const toml::table& t2) {
  auto t = t1;
  for (const auto& [k, v] : t2) {
    t.insert(k, v);
  }
  return t;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity, misc-use-internal-linkage)
int validate_subcmd(const ValidateConfig& c) {
  try {
    int return_code = 0;
    toml::table status;

    if (c.quiet) {
      // In theory nothing should write to stdout, but better to be safe than sorry
#ifdef _WIN32
      std::ignore = std::freopen("nul", "w", stdout);  // NOLINT
#else
      std::ignore = std::freopen("/dev/null", "w", stdout);  // NOLINT
#endif
    }

    const auto is_cooler = cooler::utils::is_cooler(c.uri);
    const auto is_hic = hic::utils::is_hic_file(c.uri);
    const auto is_mcool = cooler::utils::is_multires_file(c.uri, false);
    const auto is_scool = cooler::utils::is_scool_file(c.uri, false);

    if (c.include_file_path) {
      status.insert("uri", c.uri);
    }

    if (is_cooler) {
      status.insert("format", "cool");
    }
    if (is_hic) {
      status.insert("format", "hic");
    }
    if (is_mcool) {
      status.insert("format", "mcool");
    }
    if (is_scool) {
      status.insert("format", "scool");
    }

    if (!is_hic && !is_cooler && !is_mcool && !is_scool) {
      if (!c.quiet) {
        print_report(status, c.output_format);
        fmt::print(stderr, FMT_STRING("### FAILURE: \"{}\" is not in .hic or .[ms]cool format!\n"),
                   c.uri);
      }
      return 1;
    }

    if (is_hic) {
      auto res = validate_hic(c.uri, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else if (is_mcool) {
      auto res = validate_mcool(c.uri, c.validate_index, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else if (is_scool) {
      auto res = validate_scool(c.uri, c.validate_index, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else {
      auto res = validate_cooler(c.uri, c.validate_index);
      return_code = res.first;
      status = merge_tables(status, res.second);
    }

    if (!c.quiet) {
      print_report(status, c.output_format);
      if (is_hic) {
        fmt::print(stderr, FMT_STRING("### {}: \"{}\" is {}a valid .hic file."),
                   return_code == 0 ? "SUCCESS" : "FAILURE", c.uri, return_code == 0 ? "" : "not ");
      } else if (is_mcool) {
        fmt::print(stderr, FMT_STRING("### {}: \"{}\" is {}a valid .mcool file."),
                   return_code == 0 ? "SUCCESS" : "FAILURE", c.uri, return_code == 0 ? "" : "not ");
      } else if (is_scool) {
        fmt::print(stderr, FMT_STRING("### {}: \"{}\" is {}a valid .scool file."),
                   return_code == 0 ? "SUCCESS" : "FAILURE", c.uri, return_code == 0 ? "" : "not ");
      } else if (std::filesystem::exists(c.uri)) {
        fmt::print(stderr, FMT_STRING("### {}: \"{}\" is {}a valid .cool file."),
                   return_code == 0 ? "SUCCESS" : "FAILURE", c.uri, return_code == 0 ? "" : "not ");
      } else {
        fmt::print(stderr, FMT_STRING("### {}: \"{}\" {} to valid Cooler."),
                   return_code == 0 ? "SUCCESS" : "FAILURE", c.uri,
                   return_code == 0 ? "points" : "does not point");
      }
    }

    return return_code;
  } catch (...) {
    if (c.quiet) {
      return 1;
    }
    throw;
  }
}

}  // namespace hictk::tools
