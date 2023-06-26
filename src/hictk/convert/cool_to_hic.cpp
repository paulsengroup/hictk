// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>

#include <boost/process/child.hpp>
#include <boost/process/search_path.hpp>

#include "hictk/fmt.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"

namespace std {
template <>
struct default_delete<FILE> {
  void operator()(FILE* file) const { fclose(file); }
};
}  // namespace std

namespace hictk::tools {

[[nodiscard]] static std::filesystem::path find_java() {
  auto java = boost::process::search_path("java");
  if (java.empty()) {
    throw std::runtime_error("unable to find java in your PATH");
  }
  return java.string();
}

[[maybe_unused]] [[nodiscard]] static std::filesystem::path find_pigz() {
  return boost::process::search_path("pigz").string();
}

[[nodiscard]] static std::vector<std::string> generate_juicer_tools_pre_args(
    const ConvertConfig& c, const std::filesystem::path& path_to_pixels,
    const std::filesystem::path& path_to_chrom_sizes) {
  const auto heap_size = c.block_cache_size == 0 ? 8000.0 : double(c.block_cache_size) / 1.0e6;
  return {fmt::format(FMT_STRING("-Xmx{:0}M"), heap_size),
          "-jar",
          c.juicer_tools_jar.string(),
          "pre",
          "-j",
          "1",
          "--threads",
          "1",
          "-t",
          c.tmp_dir.string(),
          "-n",
          "-r",
          fmt::format(FMT_STRING("{}"), fmt::join(c.resolutions, ",")),
          path_to_pixels.string(),
          c.path_to_output.string(),
          path_to_chrom_sizes.string()};
}

[[nodiscard]] static std::vector<std::string> generate_juicer_tools_add_norm_args(
    const ConvertConfig& c, const std::filesystem::path& path_to_weights) {
  const auto heap_size = c.block_cache_size == 0 ? 8000.0 : double(c.block_cache_size) / 1.0e6;
  return {fmt::format(FMT_STRING("-Xmx{:0}M"), heap_size),
          "-jar",
          c.juicer_tools_jar.string(),
          "addNorm",
          c.path_to_output.string(),
          path_to_weights.string()};
}

static void dump_chrom_sizes(const cooler::File& clr, const std::filesystem::path& dest) {
  const std::unique_ptr<FILE> f(std::fopen(dest.c_str(), "we"));

  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), dest);
  }

  fmt::print(f.get(), FMT_STRING("{:tsv}\n"), fmt::join(clr.chromosomes(), "\n"));
}

static void dump_pixels_plain(const cooler::File& clr, const std::filesystem::path& dest) {
  const std::unique_ptr<FILE> f(std::fopen(dest.c_str(), "we"));

  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), dest);
  }
  // https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
  // <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
  for (std::uint32_t i = 0; i < clr.chromosomes().size(); ++i) {
    for (std::uint32_t j = i; j < clr.chromosomes().size(); ++j) {
      auto sel =
          clr.fetch<std::uint32_t>(clr.chromosomes().at(i).name(), clr.chromosomes().at(j).name());
      std::for_each(sel.begin(), sel.end(), [&](const Pixel<std::uint32_t>& p) {
        fmt::print(f.get(), FMT_COMPILE("0\t{}\t{}\t0\t1\t{}\t{}\t1\t{}\n"),
                   p.coords.bin1.chrom().name(), p.coords.bin1.start(),
                   p.coords.bin2.chrom().name(), p.coords.bin2.start(), p.count);
        if (!bool(f)) {  // NOLINT
          throw fmt::system_error(errno, FMT_STRING("an error occurred while pixels to file {}"),
                                  dest);
        }
      });
    }
  }
}

static void dump_pixels(const cooler::File& clr, const std::filesystem::path& dest) {
  dump_pixels_plain(clr, dest);
}

static bool dump_weights(std::uint32_t resolution, std::string_view cooler_uri,
                         const std::filesystem::path& weight_file) {
  const auto clr = cooler::File::open_read_only(cooler_uri);
  assert(clr.bin_size() == resolution);

  if (!clr.has_weights("weight")) {
    return false;
  }

  const std::unique_ptr<FILE> f(std::fopen(weight_file.c_str(), "a"));
  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), weight_file);
  }

  const auto weights = (*clr.read_weights("weight"))();
  std::ptrdiff_t i0 = 0;

  for (const auto& chrom : clr.chromosomes()) {
    // TODO add GW/INTRA/INTER prefix as appropriate
    fmt::print(f.get(), FMT_STRING("vector\tICE\t{}\t{}\tBP\n"), chrom.name(), resolution);

    const auto num_bins = (chrom.size() + resolution - 1) / resolution;
    const auto i1 = i0 + static_cast<std::ptrdiff_t>(num_bins);
    std::for_each(weights.begin() + i0, weights.begin() + i1, [&](const double w) {
      std::isnan(w) ? fmt::print(f.get(), FMT_COMPILE(".\n"))
                    : fmt::print(f.get(), FMT_COMPILE("{}\n"), 1.0 / w);
      if (!bool(f)) {  // NOLINT
        throw fmt::system_error(errno, FMT_STRING("an error occurred while weights to file {}"),
                                weight_file);
      }
    });

    i0 = i1;
  }
  return true;
}

static bool dump_weights(const ConvertConfig& c, const std::filesystem::path& weight_file) {
  bool cooler_has_weights = false;
  for (const auto& res : c.resolutions) {
    cooler_has_weights |= dump_weights(
        res, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_input.string(), res),
        weight_file);
  }

  return cooler_has_weights;
}

[[nodiscard]] static std::unique_ptr<boost::process::child> run_juicer_tools_pre(
    const ConvertConfig& c, const std::filesystem::path& chrom_sizes,
    const std::filesystem::path& pixels) {
  const auto cmd = generate_juicer_tools_pre_args(c, pixels, chrom_sizes);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

[[nodiscard]] static std::unique_ptr<boost::process::child> run_juicer_tools_add_norm(
    const ConvertConfig& c, const std::filesystem::path& path_to_weights) {
  const auto cmd = generate_juicer_tools_add_norm_args(c, path_to_weights);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

void cool_to_hic(const ConvertConfig& c) {
  static const internal::TmpDir tmpdir{};

  const auto chrom_sizes = tmpdir() / "reference.chrom.sizes";
  const auto pixels = tmpdir() / "pixels.tsv";
  const auto weights = tmpdir() / "weights.txt";

  std::unique_ptr<boost::process::child> process{};

  try {
    {
      const auto uri = c.input_format == "cool"
                           ? c.path_to_input.string()
                           : fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                         c.path_to_input.string(), c.resolutions.front());

      const auto clr = cooler::File::open_read_only(uri);
      dump_chrom_sizes(clr, chrom_sizes);
      dump_pixels(clr, pixels);
    }

    process = run_juicer_tools_pre(c, chrom_sizes, pixels);
    process->wait();
    if (process->exit_code() != 0) {
      throw std::runtime_error(fmt::format(FMT_STRING("juicer_tools pre failed with exit code {}"),
                                           process->exit_code()));
    }
    process = nullptr;

    std::filesystem::remove(chrom_sizes);
    std::filesystem::remove(pixels);

    const auto weight_file_has_data =
        c.resolutions.size() == 1
            ? dump_weights(c.resolutions.front(), c.path_to_input.string(), weights)
            : dump_weights(c, weights);

    if (weight_file_has_data) {
      process = run_juicer_tools_add_norm(c, weights);
      process->wait();
      if (process->exit_code() != 0) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("juicer_tools pre failed with exit code {}"), process->exit_code()));
      }
      process = nullptr;
    }
  } catch (const std::exception& e) {
    if (process) {
      process->terminate();
    }
    throw;
  }
}
}  // namespace hictk::tools
