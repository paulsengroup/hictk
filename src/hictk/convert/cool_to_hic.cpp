// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <boost/asio/io_context.hpp>
#include <boost/asio/write.hpp>
#include <boost/process/async_pipe.hpp>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>
#include <boost/process/search_path.hpp>

#include "hictk/fmt.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"

namespace std {
template <>
struct default_delete<FILE> {
  void operator()(FILE* file) const { std::fclose(file); }  // NOLINT
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
    const std::filesystem::path& path_to_chrom_sizes, std::size_t processes) {
  assert(processes != 0);
  const auto heap_size = c.block_cache_size == 0 ? 8000.0 : double(c.block_cache_size) / 1.0e6;
  return {fmt::format(FMT_STRING("-Xmx{:0}M"), heap_size),
          "-jar",
          c.juicer_tools_jar.string(),
          "pre",
          "-j",
          fmt::to_string(processes),
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
    const ConvertConfig& c, const std::filesystem::path& path_to_weights, std::size_t processes) {
  assert(processes != 0);
  const auto heap_size = c.block_cache_size == 0 ? 8000.0 : double(c.block_cache_size) / 1.0e6;
  return {fmt::format(FMT_STRING("-Xmx{:0}M"), heap_size),
          "-jar",
          c.juicer_tools_jar.string(),
          "addNorm",
          "-j",
          fmt::to_string(processes),
          c.path_to_output.string(),
          path_to_weights.string()};
}

static void dump_chrom_sizes(const cooler::File& clr, const std::filesystem::path& dest) {
  spdlog::info(FMT_STRING("writing chromosomes to file {}..."), dest);
  const std::unique_ptr<FILE> f(std::fopen(dest.string().c_str(), "we"));

  if (!bool(f)) {
    throw fmt::system_error(errno, FMT_STRING("cannot open file {}"), dest);
  }

  fmt::print(f.get(), FMT_STRING("{:tsv}\n"), fmt::join(clr.chromosomes(), "\n"));
  spdlog::info(FMT_STRING("DONE! Wrote {} chromosomes to file {}"), clr.chromosomes().size(), dest);
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
              spdlog::info(FMT_STRING("processing {:ucsc} {:ucsc} at {:.0f} pixels/s..."), bin1,
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
    std::size_t processes) {
  assert(compression_lvl != 0);
  assert(processes != 0);
  // clang-format off
  return std::make_unique<boost::process::child>(
      find_pigz().string(),
      fmt::format(FMT_STRING("-{}"), compression_lvl),
      "--processes", fmt::to_string(processes),
      boost::process::std_in < pipe,
      boost::process::std_out > dest.string()
  );
  // clang-format on
}

static std::size_t dump_pixels_pigz(const cooler::File& clr, const std::filesystem::path& dest,
                                    std::uint8_t compression_lvl, std::size_t processes,
                                    std::size_t update_frequency = 10'000'000) {
  assert(compression_lvl != 0);
  assert(processes > 1);

  boost::asio::io_context ioc;
  boost::process::async_pipe pipe{ioc};
  const auto pigz = run_pigz(pipe, dest, compression_lvl, processes - 1);

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
            buffer =
                fmt::format(FMT_COMPILE("0\t{}\t{}\t0\t1\t{}\t{}\t1\t{}\n"), bin1.chrom().name(),
                            bin1.start(), bin2.chrom().name(), bin2.start(), p.count);

            boost::asio::write(pipe, boost::asio::buffer(buffer.data(), buffer.size()));

            if (!pigz->running()) {
              throw std::runtime_error(fmt::format(
                  FMT_STRING("pigz returned prematurely with code {} while writing pixels to {}"),
                  pigz->exit_code(), dest));
            }
            if (++i == update_frequency) {
              const auto t1 = std::chrono::steady_clock::now();
              const auto delta =
                  static_cast<double>(
                      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
                  1000.0;
              spdlog::info(FMT_STRING("processing {:ucsc} {:ucsc} at {:.0f} pixels/s..."), bin1,
                           bin2, double(update_frequency) / delta);
              t0 = t1;
              i = 0;
            }
          });
    }
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
                        std::uint8_t compression_lvl, std::size_t processes) {
  const auto t0 = std::chrono::steady_clock::now();

  spdlog::info(FMT_STRING("writing pixels to file {}..."), dest);

  std::size_t pixels_processed{};
  if (dest.extension() == ".gz") {
    assert(compression_lvl != 0);
    pixels_processed = dump_pixels_pigz(clr, dest, compression_lvl, processes);
  } else {
    pixels_processed = dump_pixels_plain(clr, dest);
  }
  const auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()) /
      1.0e6;
  spdlog::info(FMT_STRING("wrote {} pixels across {} chromosomes to {} in {:.2f}s"),
               pixels_processed, clr.chromosomes().size(), dest, delta);
}

static bool dump_weights(std::uint32_t resolution, std::string_view cooler_uri,
                         const std::filesystem::path& weight_file) {
  spdlog::info(FMT_STRING("[{}] writing balancing weights to file {}..."), resolution, weight_file);
  const auto clr = cooler::File::open_read_only(cooler_uri);
  assert(clr.bin_size() == resolution);

  if (!clr.has_weights("weight")) {
    spdlog::warn(FMT_STRING("[{}] unable to read weights from \"{}\"..."), resolution, cooler_uri);
    return false;
  }

  const std::unique_ptr<FILE> f(std::fopen(weight_file.string().c_str(), "ae"));
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
        throw fmt::system_error(
            errno, FMT_STRING("an error occurred while writing weights to file {}"), weight_file);
      }
    });

    i0 = i1;
  }
  spdlog::info(FMT_STRING("[{}] wrote {} weights to file {}..."), resolution, weights.size(),
               weight_file);
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
    const std::filesystem::path& pixels, std::size_t processes) {
  const auto cmd = generate_juicer_tools_pre_args(c, pixels, chrom_sizes, processes);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

[[nodiscard]] static std::unique_ptr<boost::process::child> run_juicer_tools_add_norm(
    const ConvertConfig& c, const std::filesystem::path& path_to_weights, std::size_t processes) {
  const auto cmd = generate_juicer_tools_add_norm_args(c, path_to_weights, processes);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

void cool_to_hic(const ConvertConfig& c) {
  static const internal::TmpDir tmpdir{};

  const auto chrom_sizes = tmpdir() / "reference.chrom.sizes";
  const auto pixels = [&]() {
    if (c.gzip_compression_lvl == 0 || find_pigz().empty()) {
      return tmpdir() / "pixels.tsv";
    }
    return tmpdir() / "pixels.tsv.gz";
  }();
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
      dump_pixels(clr, pixels, c.gzip_compression_lvl, c.processes);
    }

    auto t1 = std::chrono::steady_clock::now();
    spdlog::info(FMT_STRING("running juicer_tools pre..."));
    process = run_juicer_tools_pre(c, chrom_sizes, pixels, c.processes);
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
    spdlog::info(FMT_STRING("DONE! Running juicer_tools pre took {:.2f}s"), delta);

    std::filesystem::remove(chrom_sizes);
    std::filesystem::remove(pixels);

    const auto weight_file_has_data =
        c.resolutions.size() == 1
            ? dump_weights(c.resolutions.front(), c.path_to_input.string(), weights)
            : dump_weights(c, weights);

    if (weight_file_has_data) {
      t1 = std::chrono::steady_clock::now();
      spdlog::info(FMT_STRING("running juicer_tools addNorm..."));
      process = run_juicer_tools_add_norm(c, weights, c.processes);
      process->wait();
      if (process->exit_code() != 0) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("juicer_tools pre failed with exit code {}"), process->exit_code()));
      }
      t2 = std::chrono::steady_clock::now();
      delta = static_cast<double>(
                  std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()) /
              1.0e6;
      spdlog::info(FMT_STRING("DONE! Running juicer_tools addNorm took {:.2f}s"), delta);
    }
  } catch (const std::exception&) {
    if (process) {
      process->terminate();
    }
    throw;
  }
}
}  // namespace hictk::tools
