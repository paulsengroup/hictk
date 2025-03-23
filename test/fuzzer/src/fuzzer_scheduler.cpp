// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <pybind11/embed.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <boost/process/v2.hpp>
#include <cassert>
#include <cstdint>
#include <exception>
#include <future>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "hictk/fuzzer/config.hpp"
#include "hictk/fuzzer/cooler.hpp"
#include "hictk/fuzzer/tools.hpp"

namespace hictk::fuzzer {

[[nodiscard]] static std::vector<std::string> config_to_cli_args(const Config& c, std::uint16_t id,
                                                                 std::uint64_t seed) {
  std::vector<std::string> args{"launch-worker",
                                c.test_uri.string(),
                                c.reference_uri.string(),
                                "--task-id",
                                fmt::to_string(id),
                                "--1d-to-2d-query-ratio",
                                fmt::to_string(c._1d_to_2d_query_ratio),
                                "--duration",
                                fmt::to_string(c.duration),
                                "--format",
                                c.query_format,
                                "--query-length-avg",
                                fmt::to_string(c.query_relative_length_avg),
                                "--query-length-std",
                                fmt::to_string(c.query_relative_length_std),
                                "--min-query-length",
                                fmt::to_string(c.min_query_length),
                                "--max-query-length",
                                fmt::to_string(c.max_query_length),
                                "--normalization",
                                c.normalization,
                                "--seed",
                                fmt::to_string(seed),
                                "--verbosity",
                                fmt::to_string(c.verbosity)};

  if (c.resolution.has_value()) {
    args.emplace_back("--resolution");
    args.emplace_back(fmt::to_string(*c.resolution));
  }

  if (c.diagonal_band_width) {
    args.emplace_back("--diagonal-band-width");
    args.emplace_back(fmt::to_string(*c.diagonal_band_width));
  }

  return args;
}

[[nodiscard]] static std::vector<std::uint64_t> generate_seeds(std::uint64_t seed,
                                                               std::size_t num_seeds) {
  assert(num_seeds != 0);
  std::vector<std::uint32_t> buffer(num_seeds * (sizeof(std::uint64_t) / sizeof(std::uint32_t)));
  std::seed_seq{seed}.generate(buffer.begin(), buffer.end());

  std::vector<std::uint64_t> seeds(num_seeds);
  // NOLINTNEXTLINE
  std::copy(buffer.begin(), buffer.end(), reinterpret_cast<std::uint32_t*>(seeds.data()));
  return seeds;
}

[[nodiscard]] static boost::process::v2::process spawn_worker_process(const Config& c,
                                                                      std::size_t id,
                                                                      std::uint64_t seed,
                                                                      boost::asio::io_context& ctx,
                                                                      std::mutex& mtx) {
  std::mt19937_64 rand_eng(seed);

  if (c.suppress_python_warnings) {
    boost::process::v2::environment::set("PYTHONWARNINGS", "ignore");
  }

  [[maybe_unused]] const auto lck = std::scoped_lock(mtx);
  for (std::size_t i = 0; i < 10; ++i) {  // NOLINT(*-avoid-magic-numbers)
    boost::process::v2::process proc(ctx, c.exec.string(),
                                     config_to_cli_args(c, static_cast<std::uint16_t>(id), seed),
                                     boost::process::v2::process_stdio{nullptr, {}, {}});

    if (proc.running()) {
      return proc;
    }
    SPDLOG_WARN(FMT_STRING("failed to spawn worker process #{} (attempt {}/10)"), id, i + 1);

    const auto sleep_time =
        std::chrono::milliseconds(std::uniform_int_distribution<std::uint64_t>{10, 500}(rand_eng));
    SPDLOG_DEBUG(
        FMT_STRING("sleeping for {}ms before attempting to launch process #{} one more time..."),
        sleep_time.count(), id);
    std::this_thread::sleep_for(sleep_time);
  }
  throw std::runtime_error(fmt::format(FMT_STRING("failed to spawn worker process #{}"), id));
}

int fuzz_subcommand(const Config& c) {
  assert(c.task_id == 0);
  [[maybe_unused]] const pybind11::scoped_interpreter guard{};
  try {
    SPDLOG_INFO(FMT_STRING("[executor] cooler version: {}"), cooler::version());

    BS::light_thread_pool tpool(c.nproc);

    assert(!c.exec.empty());
    assert(c.seed.has_value());
    std::vector<std::future<int>> futures(c.nproc);
    boost::asio::io_context ctx;
    std::mutex ctx_mtx;
    const auto seeds = generate_seeds(*c.seed, c.nproc);

    for (std::size_t i = 0; i < seeds.size(); ++i) {
      futures[i] = tpool.submit_task([&, id = i + 1, seed = seeds[i]]() {
        try {
          // NOLINTBEGIN(clang-analyzer-unix.BlockInCriticalSection)
          auto proc = spawn_worker_process(c, id, seed, ctx, ctx_mtx);
          return proc.wait();
          // NOLINTEND(clang-analyzer-unix.BlockInCriticalSection)
        } catch (const std::exception& e) {
          SPDLOG_ERROR(FMT_STRING("[{}] error occurred in worker process: {}"), id, e.what());
          return 1;
        } catch (...) {
          SPDLOG_ERROR(FMT_STRING("[{}] an unknown error occurred in worker process"), id);
          return 1;
        }
      });
    }

    int exit_code = 0;
    for (std::size_t i = 0; i < futures.size(); ++i) {
      if (const auto ec = futures[i].get(); ec != 0) {
        SPDLOG_ERROR(FMT_STRING("[{}] worker process returned exit code {}"), i + 1, ec);
        exit_code = 1;
      }
    }

    if (exit_code != 0) {
      SPDLOG_ERROR("[executor] one or more worker process returned with non-zero exit code");
    }

    return exit_code;
  } catch (const std::exception& e) {
    // wrap python exceptions while we still hold the scoped_interpreter guard
    throw std::runtime_error(e.what());
  }
}

}  // namespace hictk::fuzzer
