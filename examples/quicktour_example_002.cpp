// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cstdint>
#include <hictk/file.hpp>
#include <hictk/transformers.hpp>
#include <iostream>
#include <string>

int main() {
  const std::string path = "interactions.hic";
  const std::uint32_t resolution = 1000;

  const hictk::File f(path, resolution);

  const auto selector = f.fetch("chr1", "chr2");
  const hictk::transformers::JoinGenomicCoords jselector(
      selector.template begin<std::int32_t>(), selector.template end<std::int32_t>(), f.bins_ptr());

  for (const auto& p : jselector) {
    std::cout << p.coords.bin1.chrom().name() << "\t";
    std::cout << p.coords.bin1.start() << "\t";
    std::cout << p.coords.bin1.end() << "\t";
    std::cout << p.coords.bin2.chrom().name() << "\t";
    std::cout << p.coords.bin2.start() << "\t";
    std::cout << p.coords.bin2.end() << "\t";
    std::cout << p.count << "\n";
  }
}
