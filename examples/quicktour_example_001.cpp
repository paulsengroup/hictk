// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstdint>
#include <hictk/file.hpp>
#include <iostream>
#include <string>

int main() {
  // const std::string path = "interactions.cool";
  // const std::string path = "interactions.mcool::/resolutions/1000";
  const std::string path = "interactions.hic";
  const std::uint32_t resolution = 1000;

  const hictk::File f(path, resolution);

  const auto selector = f.fetch("chr1", "chr2");

  std::for_each(selector.template begin<std::int32_t>(), selector.template end<std::int32_t>(),
                [](const hictk::ThinPixel<std::int32_t>& p) {
                  std::cout << p.bin1_id << "\t";
                  std::cout << p.bin2_id << "\t";
                  std::cout << p.count << "\n";
                });
}
