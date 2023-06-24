// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./common.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void convert_subcmd(const ConvertConfig& c) {
  if (c.input_format == "hic") {
    hic_to_cool(c);
  } else {
    cool_to_hic(c);
  }
}
}  // namespace hictk::tools
