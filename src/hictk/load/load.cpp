// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./load.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <string>
#include <variant>

#include "hictk/pixel.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/load.hpp"

namespace hictk::tools {

[[nodiscard]] static Stats load_hic(const LoadConfig& c,  // NOLINTNEXTLINE(*-avoid-magic-numbers)
                                    std::size_t queue_capacity_bytes = 64'000'000) {
  assert(c.output_format == "hic");

  BS::light_thread_pool tpool(2);
  std::atomic<bool> early_return{false};

  const auto format = format_from_string(c.format);

  auto parser = init_pixel_parser(format, c.input_path, c.path_to_chrom_sizes, c.path_to_bin_table,
                                  c.bin_size, c.assembly, c.drop_unknown_chroms);

  const auto& bins = parser.bins();

  if (c.output_format == "hic" && bins.type() == BinTable::Type::variable) {
    throw std::runtime_error("creating a .hic file with variable bin size is not supported");
  }

  const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<float>);
  PixelQueue<float> pixel_queue{queue_capacity};

  auto producer = spawn_producer(tpool, parser, pixel_queue, c.offset, early_return,
                                 c.transpose_lower_triangular_pixels);
  auto consumer =
      spawn_consumer(tpool, c, bins, parser.assembly(), format, pixel_queue, early_return);

  try {
    producer.get();
    return consumer.get();
  } catch (...) {  // NOLINT
    tpool.wait();
    throw;
  }
}

[[nodiscard]] static Stats load_cool_float(
    const LoadConfig& c,  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::size_t queue_capacity_bytes = 64'000'000) {
  assert(c.count_as_float);
  assert(c.output_format == "cool");

  BS::light_thread_pool tpool(2);
  std::atomic<bool> early_return{false};

  const auto format = format_from_string(c.format);

  auto parser = init_pixel_parser(format, c.input_path, c.path_to_chrom_sizes, c.path_to_bin_table,
                                  c.bin_size, c.assembly, c.drop_unknown_chroms);

  const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<double>);
  PixelQueue<double> pixel_queue{queue_capacity};

  const auto& bins = parser.bins();

  auto producer = spawn_producer(tpool, parser, pixel_queue, c.offset, early_return,
                                 c.transpose_lower_triangular_pixels);
  auto consumer =
      spawn_consumer(tpool, c, bins, parser.assembly(), format, pixel_queue, early_return);

  try {
    producer.get();
    return consumer.get();
  } catch (...) {  // NOLINT
    tpool.wait();
    throw;
  }
}

[[nodiscard]] static Stats load_cool_int(
    const LoadConfig& c,  // NOLINTNEXTLINE(*-avoid-magic-numbers)
    std::size_t queue_capacity_bytes = 64'000'000) {
  assert(!c.count_as_float);
  assert(c.output_format == "cool");

  BS::light_thread_pool tpool(2);
  std::atomic<bool> early_return{false};

  const auto format = format_from_string(c.format);

  auto parser = init_pixel_parser(format, c.input_path, c.path_to_chrom_sizes, c.path_to_bin_table,
                                  c.bin_size, c.assembly, c.drop_unknown_chroms);

  const auto queue_capacity = queue_capacity_bytes / sizeof(ThinPixel<std::int32_t>);
  PixelQueue<std::int32_t> pixel_queue{queue_capacity};

  const auto& bins = parser.bins();

  auto producer = spawn_producer(tpool, parser, pixel_queue, c.offset, early_return,
                                 c.transpose_lower_triangular_pixels);
  auto consumer =
      spawn_consumer(tpool, c, bins, parser.assembly(), format, pixel_queue, early_return);

  try {
    producer.get();
    return consumer.get();
  } catch (...) {  // NOLINT
    tpool.wait();
    throw;
  }
}

int run_subcmd(const LoadConfig& c) {  // NOLINT(misc-use-internal-linkage)
  const auto t0 = std::chrono::system_clock::now();

  const auto stats = [&]() {
    if (c.count_as_float && c.output_format == "cool") {
      return load_cool_float(c);
    }
    if (!c.count_as_float && c.output_format == "cool") {
      return load_cool_int(c);
    }
    assert(c.output_format == "hic");
    return load_hic(c);
  }();

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

  std::visit(
      [&](const auto& sum) {
        SPDLOG_INFO(FMT_STRING("ingested {} interactions ({} nnz) in {}s!"), sum, stats.nnz,
                    static_cast<double>(delta) / 1.0e9);  // NOLINT(*-avoid-magic-numbers)
      },
      stats.sum);

  return 0;
}

}  // namespace hictk::tools
