// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cassert>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/tools/common.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void cool_to_hic(const ConvertConfig& c) {
  if (c.force && std::filesystem::exists(c.path_to_output)) {
    [[maybe_unused]] std::error_code ec{};
    std::filesystem::remove(c.path_to_output, ec);
  }

  const auto base_uri = c.input_format == "cool"
                            ? c.path_to_input.string()
                            : fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                          c.path_to_input.string(), c.resolutions.front());

  const cooler::File base_clr{base_uri};
  if (base_clr.bin_size() == 0) {
    throw std::runtime_error("converting cooler files with variable bin size is not supported.");
  }

  const auto chromosomes = cooler::File(base_uri).chromosomes();
  const auto& resolutions = c.resolutions;

  const internal::TmpDir tmpdir{c.tmp_dir};
  hictk::hic::internal::HiCFileWriter w(c.path_to_output.string(), chromosomes, resolutions,
                                        c.genome, c.threads, c.chunk_size, c.tmp_dir,
                                        c.compression_lvl);
  if (c.input_format == "cool") {
    w.add_pixels(base_clr.bin_size(), base_clr.begin<float>(), base_clr.end<float>());
  } else {
    assert(c.input_format == "mcool");
    const cooler::MultiResFile mclr(c.path_to_input.string());

    for (const auto& res : c.resolutions) {
      try {
        const auto clr = mclr.open(res);
        w.add_pixels(res, clr.begin<float>(), clr.end<float>());
      } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        const auto pos = msg.find("does not have interactions for resolution");
        if (pos == std::string_view::npos) {
          throw;
        }
      }
    }
  }

  w.serialize();
}
}  // namespace hictk::tools
