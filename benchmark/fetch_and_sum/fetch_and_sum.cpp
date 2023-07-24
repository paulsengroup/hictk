// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <chrono>
#include <hictk/cooler.hpp>
#include <hictk/hic.hpp>
#include <hictk/hic/utils.hpp>
#include <string_view>
#include <variant>

struct Config {
  std::string path{};
  std::string weights{"NONE"};

  std::uint32_t resolution{};
};

using namespace hictk;

[[nodiscard]] inline std::pair<std::string, std::string> parse_bedpe(std::string_view line) {
  auto parse_bed = [&]() {
    assert(!line.empty());
    const auto pos1 = line.find('\t');
    const auto pos2 = line.find('\t', pos1 + 1);
    const auto pos3 = line.find('\t', pos2 + 1);

    auto tok = std::string{line.substr(0, pos3)};
    tok[pos1] = ':';
    tok[pos2] = '-';
    line.remove_prefix(pos3 + 1);
    return tok;
  };

  return std::make_pair(parse_bed(), parse_bed());
}

template <typename PixelIt>
[[nodiscard]] static std::pair<std::size_t, double> accumulate_interactions(PixelIt first_pixel,
                                                                            PixelIt last_pixel) {
  std::size_t nnz = 0;
  const auto sum = std::accumulate(first_pixel, last_pixel, 0.0,
                                   [&](const double accumulator, const auto &pixel) {
                                     ++nnz;
                                     return accumulator + double(pixel.count);
                                   });
  return std::make_pair(nnz, sum);
}

void fetch_and_sum(const Config &c, cooler::File &&clr) {
  auto weights = clr.read_weights(c.weights);

  std::string line;
  while (std::getline(std::cin, line)) {
    const auto [range1, range2] = parse_bedpe(line);
    const auto t0 = std::chrono::system_clock::now();
    auto sel = clr.fetch(range1, range2, weights);
    const auto [nnz, sum] = accumulate_interactions(sel.begin<double>(), sel.end<double>());
    const auto t1 = std::chrono::system_clock::now();

    const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    fmt::print(FMT_STRING("{}\t{}\t{}\t{}\n"), line, nnz, sum, double(delta) / 1.0e9);
  }
}

void fetch_and_sum(const Config &c, hic::File &&hf) {
  hf.optimize_cache_size_for_random_access();
  const auto norm = balancing::Method(c.weights);

  std::string line;
  while (std::getline(std::cin, line)) {
    const auto [range1, range2] = parse_bedpe(line);
    const auto t0 = std::chrono::system_clock::now();
    auto sel = hf.fetch(range1, range2, norm);
    const auto [nnz, sum] = accumulate_interactions(sel.begin<double>(false), sel.end<double>());
    const auto t1 = std::chrono::system_clock::now();

    const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

    fmt::print(FMT_STRING("{}\t{}\t{}\t{}\n"), line, nnz, sum, double(delta) / 1.0e9);
  }
}

void fetch_and_sum(const Config &c) {
  fmt::print(FMT_STRING("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tnnz\tsum\ttime\n"));
  if (hic::utils::is_hic_file(c.path)) {
    fetch_and_sum(c, hic::File(c.path, c.resolution));
  } else {
    fetch_and_sum(c, cooler::File(c.path));
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  CLI::App cli{};
  Config config{};
  cli.add_option("file", config.path, "Path to a .cool or .hic file (Cooler URI syntax supported).")
      ->required();

  cli.add_option("--weights", config.weights,
                 "Name of the balancing weights to apply to interactions.");

  cli.add_option("--resolution", config.resolution,
                 "Matrix resolution. Ignored when input file is in Cooler format.");
  try {
    cli.parse(argc, argv);

    std::ios::sync_with_stdio(false);
    if (!config.path.empty()) {
      fetch_and_sum(config);
    }

  } catch (const CLI::ParseError &e) {
    return cli.exit(e);
  } catch (const std::exception &e) {
    assert(cli);
    fmt::print(stderr, FMT_STRING("FAILURE! {} encountered the following error: {}."), argv[0],
               e.what());
    return 1;
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! {} encountered the following error: Caught an "
                          "unhandled exception!\n"),
               argv[0]);
    return 1;
  }
  return 0;
}
