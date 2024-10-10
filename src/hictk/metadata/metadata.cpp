// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./metadata.hpp"

#include <cassert>

#include "hictk/tools/config.hpp"

namespace hictk::tools {

int metadata_subcmd(const MetadataConfig& c) {  // NOLINT(misc-use-internal-linkage)
  const auto output_format = parse_output_format(c.output_format);
  if (c.input_format == "hic") {
    return print_hic_metadata(c.uri, output_format, c.include_file_path, c.recursive);
  }
  if (c.input_format == "mcool") {
    return print_mcool_metadata(c.uri, output_format, c.include_file_path, c.recursive);
  }
  if (c.input_format == "scool") {
    return print_scool_metadata(c.uri, output_format, c.include_file_path, c.recursive);
  }
  assert(c.input_format == "cool");
  return print_cool_metadata(c.uri, output_format, c.include_file_path);
}

}  // namespace hictk::tools
