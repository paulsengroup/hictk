// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "./common.hpp"
#include "./load_cooler.hpp"
#include "./load_hic.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/compressed_reader.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

struct Header {
  std::string assembly{"unknown"};
  std::optional<Reference> chromosomes{};
};

template <Format format>
class PixelParser {
  io::CompressedReader _reader{};
  std::string _strbuff{};
  BinTable _bins{};
  std::string _assembly{};

 public:
  PixelParser() = default;
  explicit PixelParser(const std::filesystem::path& path_, std::uint32_t resolution,
                       std::string_view assembly = "unknown")
      : _reader(path_ != "-" ? io::CompressedReader{path_} : io::CompressedReader{}) {
    auto header = parse_header();
    if (!header.chromosomes) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to read chromosomes from {}"), path()));
    }

    if (header.assembly.empty() || assembly != "unknown") {
      _assembly = std::string{assembly};
    } else {
      _assembly = std::move(header.assembly);
    }

    assert(resolution != 0);
    _bins = BinTable(std::move(*header.chromosomes), resolution);
  }

  PixelParser(const std::filesystem::path& path_, BinTable bins_,
              std::string_view assembly_ = "unknown")
      : _reader(path_ != "-" ? io::CompressedReader{path_} : io::CompressedReader{}),
        _bins(std::move(bins_)),
        _assembly(std::string{assembly_}) {
    try {
      std::ignore = parse_header();
    } catch (const std::exception& e) {
      SPDLOG_WARN(FMT_STRING("encountered an error while parsing the file header: {}"), e.what());
    }
  }

  [[nodiscard]] const std::filesystem::path& path() const noexcept {
    static const std::filesystem::path stdin_{"stdin"};
    if (!!_reader) {
      return stdin_;
    }
    return _reader.path();
  }

  [[nodiscard]] std::string_view assembly() const noexcept { return _assembly; }

  [[nodiscard]] const BinTable& bins() const noexcept { return _bins; }

  template <typename N>
  [[nodiscard]] ThinPixel<N> next_pixel(std::int64_t offset) {
    ThinPixel<N> p{};
    std::ignore = next_pixel(p, offset);
    return p;
  }

  template <typename N>
  [[nodiscard]] bool next_pixel(ThinPixel<N>& buff, std::int64_t offset) {
    if (_strbuff.empty()) {  // EOF
      buff.bin1_id = ThinPixel<N>::null_id;
      buff.bin2_id = ThinPixel<N>::null_id;
      buff.count = 0;
      return false;
    }

    switch (format) {
      case Format::COO:
        buff = ThinPixel<N>::from_coo(bins(), _strbuff, offset);
        break;
      case Format::BG2:
        buff = Pixel<N>::from_bg2(bins(), _strbuff, offset).to_thin();
        break;
      case Format::VP:
        buff = Pixel<N>::from_validpair(bins(), _strbuff, offset).to_thin();
        break;
      case Format::_4DN:
        buff = Pixel<N>::from_4dn_pairs(bins(), _strbuff, offset).to_thin();
        break;
    }
    std::ignore = getline();

    return true;
  }

 private:
  [[nodiscard]] bool getline(char sep = '\n') {
    const auto ok =
        !_reader ? !!std::getline(std::cin, _strbuff, sep) : _reader.getline(_strbuff, sep);
    if (!ok) {
      return ok;
    }

    if (_strbuff.empty()) {
      return getline(sep);
    }

    if (_strbuff.back() == '\r') {
      _strbuff.pop_back();
    }

    return ok;
  }

  [[nodiscard]] static std::string_view remove_prefix(std::string_view s, std::string_view prefix,
                                                      bool strip_whitespaces = true) noexcept {
    if (internal::starts_with(s, prefix)) {
      s.remove_prefix(prefix.size());
    }

    while (strip_whitespaces && !s.empty()) {
      if (std::isspace(s.front())) {
        s.remove_prefix(1);
      }
    }

    return s;  // NOLINT
  }

  [[nodiscard]] static std::pair<std::string, std::uint32_t> parse_chromsize(
      std::string_view line) {
    assert(internal::starts_with(line, "#chromsize:"));
    line = remove_prefix(line, "#chromsize:");
    auto it = std::find_if(line.begin(), line.end(), [](const char c) { return std::isspace(c); });
    if (it == line.begin() || it == line.end()) {
      throw std::runtime_error(fmt::format(FMT_STRING("malformed chromsize entry \"{}\"."), line));
    }
    std::string chrom_name{line.begin(), it};

    it = std::find_if(line.begin(), line.end(), [](const char c) { return !std::isspace(c); });
    const auto strlen = static_cast<std::size_t>(std::distance(it, line.end()));

    std::string_view chrom_size{it, strlen};
    if (chrom_size.empty()) {
      throw std::runtime_error(fmt::format(FMT_STRING("malformed chromsize entry \"{}\"."), line));
    }

    try {
      return std::make_pair(std::move(chrom_name),
                            internal::parse_numeric_or_throw<std::uint32_t>(chrom_size));
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("malformed chromsize entry \"{}\": {}"), line, e.what()));
    }
  }

  [[nodiscard]] Header parse_header() {
    if constexpr (format != Format::_4DN) {
      return {};
    }

    std::optional<std::string> assembly{};
    std::vector<std::string> chrom_names{};
    std::vector<std::uint32_t> chrom_sizes{};
    for (std::size_t i = 0; getline(); ++i) {
      try {
        if (i == 0 && _strbuff != "## pairs format v1.0") {
          throw std::runtime_error(
              "invalid header: first line in input file does not start with \"## pairs format "
              "v1.0\"");
        }
        if (_strbuff.empty()) {
          continue;
        }

        if (_strbuff.front() != '#') {
          break;
        }

        if (internal::starts_with(_strbuff, "#genome_assembly:")) {
          if (assembly.has_value()) {
            throw std::runtime_error(
                "found duplicate entry for \"genome_assembly\" in file header.");
          }
          assembly = std::string{remove_prefix(_strbuff, "#genome_assembly:")};
          continue;
        }

        if (!internal::starts_with(_strbuff, "#chromsize:")) {
          continue;
        }
        auto [chrom_name, chrom_size] = parse_chromsize(_strbuff);

        chrom_names.emplace_back(std::move(chrom_name));
        chrom_sizes.push_back(chrom_size);

      } catch (const std::exception& e) {
        throw std::runtime_error(fmt::format(FMT_STRING("failed to parse line {} from {}: {}"),
                                             i + 1, path(), e.what()));
      }
    }

    if (!assembly.has_value()) {
      assembly = "unknown";
    }

    return {*assembly, Reference{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()}};
  }
};

using PixelParserVar = std::variant<PixelParser<Format::_4DN>, PixelParser<Format::BG2>,
                                    PixelParser<Format::COO>, PixelParser<Format::VP>>;

template <typename N, Format format>
static void parse_pixels(PixelParser<format>& parser, std::int64_t offset, PixelQueue<N>& queue,
                         std::atomic<bool>& early_return) {
  ThinPixel<N> buffer{};
  while (!early_return && parser.next_pixel(buffer, offset)) {
    assert(buffer.bin1_id != ThinPixel<N>::null_id);
    assert(buffer.bin2_id != ThinPixel<N>::null_id);
    assert(buffer.count != 0);
    while (!queue.try_enqueue(buffer) && !early_return) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }

  if (!early_return) {
    buffer.bin1_id = ThinPixel<N>::null_id;
    buffer.bin2_id = ThinPixel<N>::null_id;
    buffer.count = 0;
    while (!queue.try_enqueue(buffer) && !early_return) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }
}

[[nodiscard]] static BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                             std::uint32_t bin_size) {
  auto chroms = Reference::from_chrom_sizes(path_to_chrom_sizes);
  return {chroms, bin_size};
}

[[nodiscard]] static BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                             const std::filesystem::path& path_to_bin_table) {
  auto chroms = Reference::from_chrom_sizes(path_to_chrom_sizes);

  std::ifstream ifs{};
  ifs.exceptions(std::ios::badbit);
  ifs.open(path_to_bin_table);

  std::vector<std::uint32_t> start_pos{};
  std::vector<std::uint32_t> end_pos{};

  std::string line{};
  GenomicInterval record{};
  bool fixed_bin_size = true;
  std::uint32_t bin_size = 0;
  while (std::getline(ifs, line)) {
    record = GenomicInterval::parse_bed(chroms, line);
    if (bin_size == 0) {
      bin_size = record.size();
    }

    fixed_bin_size &= record.size() == bin_size || record.chrom().size() == record.end();

    start_pos.push_back(record.start());
    end_pos.push_back(record.end());
  }

  if (fixed_bin_size) {
    SPDLOG_INFO(FMT_STRING("detected bin table with uniform bin size."));
    return {chroms, bin_size};
  }

  SPDLOG_INFO(FMT_STRING("detected bin table with variable bin size."));
  return {chroms, start_pos, end_pos};
}

[[nodiscard]] static BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                             const std::filesystem::path& path_to_bin_table,
                                             std::uint32_t bin_size) {
  if (!path_to_bin_table.empty()) {
    return init_bin_table(path_to_chrom_sizes, path_to_bin_table);
  }
  return init_bin_table(path_to_chrom_sizes, bin_size);
}

[[nodiscard]] static PixelParserVar init_pixel_parser(
    const Format format, const std::filesystem::path& path_to_interactions,
    const std::filesystem::path& path_to_chrom_sizes, const std::filesystem::path& path_to_bins,
    std::uint32_t resolution, std::string_view assembly) {
  auto bins = path_to_chrom_sizes.empty()
                  ? BinTable{}
                  : init_bin_table(path_to_chrom_sizes, path_to_bins, resolution);
  switch (format) {
    case Format::_4DN: {
      if (path_to_chrom_sizes.empty()) {
        assert(resolution != 0);
        return PixelParser<Format::_4DN>{path_to_interactions, resolution, assembly};
      }
      return PixelParser<Format::_4DN>{path_to_interactions, std::move(bins), assembly};
    }
    case Format::BG2: {
      assert(!path_to_chrom_sizes.empty());
      return PixelParser<Format::BG2>{path_to_interactions, std::move(bins), assembly};
    }
    case Format::COO: {
      assert(!path_to_chrom_sizes.empty());
      return PixelParser<Format::COO>{path_to_interactions, std::move(bins), assembly};
    }
    case Format::VP: {
      assert(!path_to_chrom_sizes.empty());
      return PixelParser<Format::VP>{path_to_interactions, std::move(bins), assembly};
    }
  }
}

template <std::size_t queue_capacity_bytes = 64'000'000>
[[nodiscard]] static PixelQueueVar init_pixel_queue(std::string_view output_format,
                                                    bool count_as_float) {
  if (count_as_float && output_format == "cool") {
    const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<double>);
    return {PixelQueue<double>{queue_capacity}};
  }
  if (count_as_float || output_format == "hic") {
    const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<float>);
    return {PixelQueue<float>{queue_capacity}};
  }

  assert(!count_as_float && output_format == "cool");
  const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<std::int32_t>);
  return {PixelQueue<std::int32_t>{queue_capacity}};
}

static Stats ingest_pixels_hic(const LoadConfig& c, const BinTable& bins,
                               const std::string& assembly, PixelQueue<float>& queue,
                               const std::atomic<bool>& early_return) {
  assert(c.output_format == "hic");
  const internal::TmpDir tmpdir{c.tmp_dir, true};
  return ingest_pixels_hic(queue, early_return, c.output_path, tmpdir(), bins.chromosomes(),
                           bins.resolution(), assembly, c.skip_all_vs_all_matrix, c.threads,
                           c.batch_size, c.compression_lvl, c.force);
}

template <typename N>
static Stats ingest_pixels_cooler(const LoadConfig& c, const BinTable& bins,
                                  std::string_view assembly, PixelQueue<N>& queue,
                                  const std::atomic<bool>& early_return) {
  assert(c.output_format == "cool");
  const internal::TmpDir tmpdir{c.tmp_dir, true};
  const auto tmp_cooler_path =
      (tmpdir() / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

  return c.assume_sorted ? ingest_pixels_sorted_cooler(queue, early_return, c.output_path, bins,
                                                       assembly, c.batch_size, c.compression_lvl,
                                                       c.force, c.validate_pixels)
                         : ingest_pixels_unsorted_cooler(
                               queue, early_return, c.output_path, tmp_cooler_path, bins, assembly,
                               c.batch_size, c.compression_lvl, c.force, c.validate_pixels);
}

template <typename N>
static Stats ingest_pairs_cooler(const LoadConfig& c, const BinTable& bins,
                                 std::string_view assembly, PixelQueue<N>& queue,
                                 const std::atomic<bool>& early_return) {
  const internal::TmpDir tmpdir{c.tmp_dir, true};
  const auto tmp_cooler_path =
      (tmpdir() / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

  return ingest_pairs_cooler(queue, early_return, c.output_path, tmp_cooler_path, bins, assembly,
                             c.batch_size, c.compression_lvl, c.force, c.validate_pixels);
}

static Stats ingest_pairs_hic(const LoadConfig& c, const BinTable& bins,
                              const std::string& assembly, PixelQueue<float>& queue,
                              const std::atomic<bool>& early_return) {
  const internal::TmpDir tmpdir{c.tmp_dir, true};
  return ingest_pairs_hic(queue, early_return, c.output_path, tmpdir(), bins.chromosomes(),
                          bins.resolution(), assembly, c.skip_all_vs_all_matrix, c.threads,
                          c.batch_size, c.compression_lvl, c.force);
}

template <typename N>
[[nodiscard]] static Stats ingest_pixels(const LoadConfig& c, const BinTable& bins,
                                         std::string_view assembly, PixelQueue<N>& queue,
                                         const std::atomic<bool>& early_return) {
  if (c.output_format == "hic") {
    if constexpr (std::is_same_v<N, float>) {
      return ingest_pixels_hic(c, bins, std::string{assembly}, queue, early_return);
    } else {
      throw std::logic_error(
          "ingest_pixels() was called with a pixel count type different than float. This is not "
          "supported when output format is .hic!");
    }
  }

  if constexpr (!std::is_same_v<N, std::int32_t> && !std::is_same_v<N, double>) {
    throw std::logic_error(
        "ingest_pixels() was called with a pixel count type different than std::int32_t and "
        "double. This is not supported when the output format is .cool!");
  } else {
    return ingest_pixels_cooler(c, bins, assembly, queue, early_return);
  }
}

template <typename N>
[[nodiscard]] static Stats ingest_pairs(const LoadConfig& c, const BinTable& bins,
                                        std::string_view assembly, PixelQueue<N>& queue,
                                        const std::atomic<bool>& early_return) {
  if (c.output_format == "hic") {
    if constexpr (std::is_same_v<N, float>) {
      return ingest_pairs_hic(c, bins, std::string{assembly}, queue, early_return);
    } else {
      throw std::logic_error(
          "ingest_pairs() was called with a pixel count type different than float. This is not "
          "supported when output format is .hic!");
    }
  }

  if constexpr (!std::is_same_v<N, std::int32_t> && !std::is_same_v<N, double>) {
    throw std::logic_error(
        "ingest_pixels() was called with a pixel count type different than std::int32_t and "
        "double. This is not supported when the output format is .cool!");
  } else {
    return ingest_pairs_cooler(c, bins, assembly, queue, early_return);
  }
}

template <typename N, Format format>
[[nodiscard]] static std::future<void> spawn_producer(BS::thread_pool& tpool,
                                                      PixelParser<format>& parser,
                                                      PixelQueue<N>& queue, std::int64_t offset,
                                                      std::atomic<bool>& early_return) {
  return tpool.submit_task([&]() {
    try {
      return parse_pixels(parser, offset, queue, early_return);
    } catch (...) {
      SPDLOG_WARN(
          FMT_STRING("exception caught in thread parsing interactions: returning immediately!"));
      early_return = true;
      throw;
    }
  });
}

template <typename N>
[[nodiscard]] static std::future<Stats> spawn_consumer(BS::thread_pool& tpool, const LoadConfig& c,
                                                       const BinTable& bins,
                                                       std::string_view assembly, Format format,
                                                       PixelQueue<N>& queue,
                                                       std::atomic<bool>& early_return) {
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;

  return tpool.submit_task([&]() -> Stats {
    try {
      return pixel_has_count ? ingest_pixels(c, bins, assembly, queue, early_return)
                             : ingest_pairs(c, bins, assembly, queue, early_return);
    } catch (...) {
      SPDLOG_WARN(FMT_STRING("exception caught in thread writing interactions to file "
                             "\"{}\": returning immediately!"),
                  c.output_path);
      early_return = true;
      throw;
    }
  });
}

int load_subcmd(const LoadConfig& c) {
  BS::thread_pool tpool(2);
  std::atomic<bool> early_return{false};

  const auto t0 = std::chrono::system_clock::now();

  const auto format = format_from_string(c.format);

  auto parser_var = init_pixel_parser(format, c.input_path, c.path_to_chrom_sizes,
                                      c.path_to_bin_table, c.bin_size, c.assembly);

  auto pixel_queue_var = init_pixel_queue(c.output_format, c.count_as_float);
  const auto stats = std::visit(
      [&](auto& parser) -> Stats {
        const auto& bins = parser.bins();

        if (c.output_format == "hic" && bins.type() == BinTable::Type::variable) {
          throw std::runtime_error("creating a .hic file with variable bin size is not supported");
        }

        return std::visit(
            [&](auto& queue) -> Stats {
              auto producer = spawn_producer(tpool, parser, queue, c.offset, early_return);
              auto consumer =
                  spawn_consumer(tpool, c, bins, parser.assembly(), format, queue, early_return);

              producer.get();
              return consumer.get();
            },
            pixel_queue_var);
      },
      parser_var);

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

  std::visit(
      [&](const auto& sum) {
        SPDLOG_INFO(FMT_STRING("ingested {} interactions ({} nnz) in {}s!"), sum, stats.nnz,
                    static_cast<double>(delta) / 1.0e9);
      },
      stats.sum);

  return 0;
}

}  // namespace hictk::tools
