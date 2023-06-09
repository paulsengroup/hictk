// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static int validate_hic(const std::string& path, bool quiet) {
  if (!quiet) {
    fmt::print(FMT_STRING("### SUCCESS: \"{}\" is a valid .hic file.\n"), path);
  }
  return 0;
}

static int validate_cooler(std::string_view path, bool validate_index, bool quiet) {
  auto status = cooler::utils::is_cooler(path);
  auto index_ok = true;
  if (validate_index) {
    index_ok = cooler::utils::index_is_valid(path, quiet);
  }

  const auto cooler_is_valid = !!status && index_ok;

  if (!quiet) {
    fmt::print(FMT_STRING("{}\n"), status);
    if (validate_index) {
      fmt::print(FMT_STRING("index_is_valid={}\n"), index_ok);
    } else {
      fmt::print(FMT_STRING("index_is_valid=not_checked\n"));
    }

    fmt::print(FMT_STRING("### {}: \"{}\" {} a valid Cooler.\n"),
               cooler_is_valid ? "SUCCESS" : "FAILURE", path, cooler_is_valid ? "is" : "is not");
  }
  return static_cast<int>(!cooler_is_valid);
}

int validate_subcmd(const ValidateConfig& c) {
  const auto is_hic = hic::utils::is_hic_file(c.uri);
  const auto is_cooler = cooler::utils::is_cooler(c.uri);

  if (!is_hic && !is_cooler) {
    if (!c.quiet) {
      fmt::print(FMT_STRING("### FAILURE: \"{}\" is not in .hic or .cool format!\n"), c.uri);
    }
    return 1;
  }

  if (is_hic) {
    return validate_hic(c.uri, c.quiet);
  }

  return validate_cooler(c.uri, c.validate_index, c.quiet);
}

}  // namespace hictk::tools
