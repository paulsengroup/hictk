// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <boost/asio/buffer.hpp>
#include <boost/asio/io_context.hpp>
#include <boost/asio/write.hpp>
#include <boost/process/async_pipe.hpp>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>
#include <boost/process/search_path.hpp>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/juicer_tools.hpp"

namespace hictk::tools {

[[maybe_unused]] [[nodiscard]] static std::filesystem::path find_pigz() {
  return boost::process::search_path("pigz").string();
}

static void dump_chrom_sizes(const cooler::File& clr, const std::filesystem::path& dest) {
  SPDLOG_INFO(FMT_STRING("writing chromosomes to file {}..."), dest);
  const std::unique_ptr<FILE> f(std::fopen(dest.string().c_str(), "we"));

  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), dest);
  }

  fmt::print(f.get(), FMT_STRING("{:tsv}\n"), fmt::join(clr.chromosomes(), "\n"));
  SPDLOG_INFO(FMT_STRING("DONE! Wrote {} chromosomes to file {}"), clr.chromosomes().size(), dest);
}

static std::size_t dump_pixels_plain(const cooler::File& clr, const std::filesystem::path& dest,
                                     std::size_t update_frequency = 10'000'000) {
  const std::unique_ptr<FILE> f(std::fopen(dest.string().c_str(), "we"));

  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), dest);
  }

  std::size_t i = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (std::uint32_t chrom1_id = 0; chrom1_id < clr.chromosomes().size(); ++chrom1_id) {
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < clr.chromosomes().size(); ++chrom2_id) {
      auto sel =
          clr.fetch(clr.chromosomes().at(chrom1_id).name(), clr.chromosomes().at(chrom2_id).name());
      std::for_each(
          sel.begin<std::int32_t>(), sel.end<std::int32_t>(),
          [&](const ThinPixel<std::int32_t>& p) {
            const auto bin1 = clr.bins().at(p.bin1_id);
            const auto bin2 = clr.bins().at(p.bin2_id);
            // https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
            // <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
            fmt::print(f.get(), FMT_COMPILE("0\t{}\t{}\t0\t1\t{}\t{}\t1\t{}\n"),
                       bin1.chrom().name(), bin1.start(), bin2.chrom().name(), bin2.start(),
                       p.count);
            if (!bool(f)) {  // NOLINT
              throw fmt::system_error(
                  errno, FMT_STRING("an error occurred while pixels to file {}"), dest);
            }

            if (++i == update_frequency) {
              const auto t1 = std::chrono::steady_clock::now();
              const auto delta =
                  static_cast<double>(
                      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
                  1000.0;
              SPDLOG_INFO(FMT_STRING("processing {:ucsc} {:ucsc} at {:.0f} pixels/s..."), bin1,
                          bin2, double(update_frequency) / delta);
              t0 = t1;
              i = 0;
            }
          });
    }
  }

  assert(clr.attributes().nnz);  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return static_cast<std::size_t>(*clr.attributes().nnz);
}

template <typename Pipe>
[[nodiscard]] static std::unique_ptr<boost::process::child> run_pigz(
    Pipe& pipe, const std::filesystem::path& dest, std::uint8_t compression_lvl,
    std::size_t threads) {
  assert(compression_lvl != 0);
  assert(threads != 0);
  // clang-format off
  return std::make_unique<boost::process::child>(
      find_pigz().string(),
      fmt::format(FMT_STRING("-{}"), compression_lvl),
      "--processes", fmt::to_string(threads),
      boost::process::std_in < pipe,
      boost::process::std_out > dest.string()
  );
  // clang-format on
}

static std::size_t dump_pixels_pigz(const cooler::File& clr, const std::filesystem::path& dest,
                                    std::uint8_t compression_lvl, std::size_t threads,
                                    std::size_t update_frequency = 10'000'000) {
  assert(compression_lvl != 0);
  assert(threads > 1);

  boost::asio::io_context ioc;
  boost::process::async_pipe pipe{ioc};
  const auto pigz = run_pigz(pipe, dest, compression_lvl, threads - 1);

  auto t0 = std::chrono::steady_clock::now();
  std::string buffer;
  std::size_t i = 0;
  for (std::uint32_t chrom1_id = 0; chrom1_id < clr.chromosomes().size(); ++chrom1_id) {
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < clr.chromosomes().size(); ++chrom2_id) {
      auto sel =
          clr.fetch(clr.chromosomes().at(chrom1_id).name(), clr.chromosomes().at(chrom2_id).name());
      std::for_each(
          sel.begin<std::int32_t>(), sel.end<std::int32_t>(),
          [&](const ThinPixel<std::int32_t>& p) {
            const auto bin1 = clr.bins().at(p.bin1_id);
            const auto bin2 = clr.bins().at(p.bin2_id);
            // https://github.com/aidenlab/juicer/wiki/Pre#short-with-score-format
            // <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <score>
            buffer +=
                fmt::format(FMT_COMPILE("0\t{}\t{}\t0\t1\t{}\t{}\t1\t{}\n"), bin1.chrom().name(),
                            bin1.start(), bin2.chrom().name(), bin2.start(), p.count);

            if (buffer.size() > 65'000) {
              if (!pigz->running()) {
                throw std::runtime_error(fmt::format(
                    FMT_STRING("pigz returned prematurely with code {} while writing pixels to {}"),
                    pigz->exit_code(), dest));
              }
              boost::asio::write(pipe, boost::asio::buffer(buffer.data(), buffer.size()));
              buffer.clear();
            }

            if (++i == update_frequency) {
              const auto t1 = std::chrono::steady_clock::now();
              const auto delta =
                  static_cast<double>(
                      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
                  1000.0;
              SPDLOG_INFO(FMT_STRING("processing {:ucsc} {:ucsc} at {:.0f} pixels/s..."), bin1,
                          bin2, double(update_frequency) / delta);
              t0 = t1;
              i = 0;
            }
          });
    }
  }

  if (!pigz->running()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("pigz returned prematurely with code {} while writing pixels to {}"),
                    pigz->exit_code(), dest));
  }
  if (!buffer.empty()) {
    boost::asio::write(pipe, boost::asio::buffer(buffer.data(), buffer.size()));
  }

  pipe.close();
  ioc.run();
  pigz->wait();
  if (pigz->exit_code() != 0) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("pigz failed with exit code {}"), pigz->exit_code()));
  }

  assert(clr.attributes().nnz);  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  return static_cast<std::size_t>(*clr.attributes().nnz);
}

static void dump_pixels(const cooler::File& clr, const std::filesystem::path& dest,
                        std::uint8_t compression_lvl, std::size_t threads) {
  const auto t0 = std::chrono::steady_clock::now();

  SPDLOG_INFO(FMT_STRING("writing pixels to file {}..."), dest);

  std::size_t pixels_processed{};
  if (dest.extension() == ".gz") {
    assert(compression_lvl != 0);
    pixels_processed = dump_pixels_pigz(clr, dest, compression_lvl, threads);
  } else {
    pixels_processed = dump_pixels_plain(clr, dest);
  }
  const auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
      1.0e6;
  SPDLOG_INFO(FMT_STRING("wrote {} pixels across {} chromosomes to {} in {:.2f}s"),
              pixels_processed, clr.chromosomes().size(), dest, delta);
}

[[nodiscard]] static std::shared_ptr<const balancing::Weights> try_read_weights(
    const cooler::File& clr, const balancing::Method& method) {
  try {
    return clr.read_weights(method);
  } catch (const std::exception& e) {
    return clr.read_weights(method, balancing::Weights::Type::DIVISIVE);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static bool dump_weights(std::uint32_t resolution, std::string_view cooler_uri,
                         const std::filesystem::path& weight_file,
                         std::vector<balancing::Method> normalizations, bool fail_if_norm_missing) {
  if (normalizations.size() == 1 && normalizations.front() == balancing::Method::NONE()) {
    return false;
  }

  if (normalizations.empty()) {
    normalizations = cooler::File(cooler_uri).avail_normalizations();
    if (normalizations.empty()) {
      return false;
    }
  }

  SPDLOG_INFO(FMT_STRING("[{}] writing balancing weights to file {}..."), resolution, weight_file);
  const cooler::File clr(cooler_uri);
  assert(clr.bin_size() == resolution);

  const std::unique_ptr<FILE> f(std::fopen(weight_file.string().c_str(), "ae"));
  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), weight_file);
  }

  for (const auto& norm : normalizations) {
    if (!clr.has_normalization(norm) && !fail_if_norm_missing) {
      SPDLOG_WARN(FMT_STRING("[{}] unable to read weights from \"{}\"..."), resolution, cooler_uri);
      continue;
    }

    const auto weights = try_read_weights(clr, norm);
    const auto weight_name = norm == "weight" ? "ICE" : norm.to_string();
    const auto weight_is_divisive = weights->type() == balancing::Weights::Type::INFER ||
                                    weights->type() == balancing::Weights::Type::UNKNOWN ||
                                    weights->type() == balancing::Weights::Type::DIVISIVE;
    auto weight_vector = (*weights)();
    if (weight_is_divisive) {
      std::transform(weight_vector.begin(), weight_vector.end(), weight_vector.begin(),
                     [](const double w) { return 1.0 / w; });
    }

    std::ptrdiff_t i0 = 0;
    for (const auto& chrom : clr.chromosomes()) {
      // TODO add GW/INTRA/INTER prefix as appropriate
      fmt::print(f.get(), FMT_STRING("vector\t{}\t{}\t{}\tBP\n"), weight_name, chrom.name(),
                 resolution);

      const auto num_bins = (chrom.size() + resolution - 1) / resolution;
      const auto i1 = i0 + static_cast<std::ptrdiff_t>(num_bins);
      std::for_each(weight_vector.begin() + i0, weight_vector.begin() + i1, [&](double w) {
        !std::isfinite(w) ? fmt::print(f.get(), FMT_COMPILE(".\n"))
                          : fmt::print(f.get(), FMT_COMPILE("{}\n"), w);
        if (!bool(f)) {  // NOLINT
          throw fmt::system_error(
              errno, FMT_STRING("an error occurred while writing weights to file {}"), weight_file);
        }
      });

      i0 = i1;
    }
    SPDLOG_INFO(FMT_STRING("[{}] wrote \"{}\" weights to file {}..."), resolution, weight_name,
                weight_file);
  }

  return std::ftell(f.get()) != 0;
}

static bool dump_weights(const ConvertConfig& c, const std::filesystem::path& weight_file) {
  bool cooler_has_weights = false;
  for (const auto& res : c.resolutions) {
    cooler_has_weights |= dump_weights(
        res, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.path_to_input.string(), res),
        weight_file, c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable);
  }

  return cooler_has_weights;
}

void cool_to_hic(const ConvertConfig& c) {
  std::ignore = find_java();

  const internal::TmpDir tmpdir{c.tmp_dir};

  const auto chrom_sizes = tmpdir() / "reference.chrom.sizes";
  const auto pixels = [&]() {
    if (c.gzip_compression_lvl == 0 || find_pigz().empty()) {
      return tmpdir() / "pixels.tsv";
    }
    return tmpdir() / "pixels.tsv.gz";
  }();
  const auto weights = tmpdir() / "weights.txt";

  if (c.force && std::filesystem::exists(c.path_to_output)) {
    [[maybe_unused]] std::error_code ec{};
    std::filesystem::remove(c.path_to_output, ec);
  }

  std::unique_ptr<boost::process::child> process{};

  try {
    {
      const auto uri = c.input_format == "cool"
                           ? c.path_to_input.string()
                           : fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                         c.path_to_input.string(), c.resolutions.front());

      const cooler::File clr(uri);
      dump_chrom_sizes(clr, chrom_sizes);
      dump_pixels(clr, pixels, c.gzip_compression_lvl, c.threads);
    }

    auto t1 = std::chrono::steady_clock::now();
    SPDLOG_INFO(FMT_STRING("running juicer_tools pre..."));
    process = run_juicer_tools_pre(c, chrom_sizes, pixels, c.threads);
    process->wait();
    if (process->exit_code() != 0) {
      throw std::runtime_error(fmt::format(FMT_STRING("juicer_tools pre failed with exit code {}"),
                                           process->exit_code()));
    }
    process = nullptr;
    auto t2 = std::chrono::steady_clock::now();
    auto delta = static_cast<double>(
                     std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) /
                 1.0e6;
    SPDLOG_INFO(FMT_STRING("DONE! Running juicer_tools pre took {:.2f}s"), delta);

    std::filesystem::remove(chrom_sizes);
    std::filesystem::remove(pixels);

    const auto weight_file_has_data =
        c.input_format == "cool"
            ? dump_weights(c.resolutions.front(), c.path_to_input.string(), weights,
                           c.normalization_methods, c.fail_if_normalization_method_is_not_avaliable)
            : dump_weights(c, weights);

    if (weight_file_has_data) {
      t1 = std::chrono::steady_clock::now();
      SPDLOG_INFO(FMT_STRING("running juicer_tools addNorm..."));
      process = run_juicer_tools_add_norm(c.juicer_tools_jar, weights, c.path_to_output,
                                          c.juicer_tools_xmx);
      process->wait();
      if (process->exit_code() != 0) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("juicer_tools pre failed with exit code {}"), process->exit_code()));
      }
      t2 = std::chrono::steady_clock::now();
      delta = static_cast<double>(
                  std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) /
              1.0e6;
      SPDLOG_INFO(FMT_STRING("DONE! Running juicer_tools addNorm took {:.2f}s"), delta);
    }
  } catch (const std::exception&) {
    if (process) {
      process->terminate();
    }
    throw;
  }
}
}  // namespace hictk::tools
