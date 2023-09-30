// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstdint>
#include <hictk/file.hpp>
#include <iostream>
#include <string>
#include <variant>

int main() {
  const std::string path = "interactions.hic";
  const std::uint32_t resolution = 1000;

  const hictk::File f(path, resolution);

  const auto selector = f.fetch("chr1", "chr2");

  // std::visit applies the lambda funciton provided as first argument
  // to the variant returned by selector.get().
  // In this way, the type held by the std::variant is checked once
  // and the underlying PixelSelector and iterators are used for all operations
  std::visit(
      [&](const auto& sel) {
        std::for_each(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(),
                      [](const hictk::ThinPixel<std::int32_t>& p) {
                        std::cout << p.bin1_id << "\t";
                        std::cout << p.bin2_id << "\t";
                        std::cout << p.count << "\n";
                      });
      },
      selector.get());
}
