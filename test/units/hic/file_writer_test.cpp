// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef SPDLOG_ACTIVE_LEVEL
#undef SPDLOG_ACTIVE_LEVEL
#endif

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE

#include "hictk/hic/file_writer.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>

#include "hictk/chromosome.hpp"
#include "hictk/hic.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"
#include "tmpdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_writer {

using namespace hictk::hic::internal;

TEST_CASE("HiC: HiCInteractionToBlockMapper::BlockMapper", "[hic][v9][short]") {
  SECTION("intra") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3R:0-50000
    {
      const HiCInteractionToBlockMapper::BlockMapperIntra mapper(803, 4);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 100) == 0);
      CHECK(mapper(802, 802) == 0);
      CHECK(mapper(803, 803) == 1);
      CHECK(mapper(1038, 2137) == 1);
      CHECK(mapper(235, 1376) == 5);
      CHECK(mapper(8, 3203) == 5);
    }
  }

  SECTION("inter") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3L:0-50000 --range2
    // chr3R:0-10000000
    {
      const HiCInteractionToBlockMapper::BlockMapperInter mapper(803, 4);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 100) == 0);
      CHECK(mapper(802, 802) == 0);
      CHECK(mapper(7, 803) == 4);
      CHECK(mapper(795, 1605) == 4);
    }

    {
      const HiCInteractionToBlockMapper::BlockMapperInter mapper(101, 1);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 99) == 0);
      CHECK(mapper(99, 99) == 0);
    }
  }
}

TEST_CASE("HiC: HiCInteractionToBlockMapper", "[hic][v9][short]") {
  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_block_partitioner.bin").string();
  const std::uint32_t resolution = 25'000;

  const hic::File f1(path1, resolution);
  const auto sel1 = f1.fetch("chr2L");
  const auto sel2 = f1.fetch("chr2L", "chr2R");

  const std::vector<ThinPixel<float>> pixels1(sel1.begin<float>(), sel1.end<float>());
  const std::vector<ThinPixel<float>> pixels2(sel2.begin<float>(), sel2.end<float>());

  HiCInteractionToBlockMapper partitioner(path2, f1.bins_ptr(), 50'000, 3);

  partitioner.append_pixels(pixels1.begin(), pixels1.end());
  partitioner.append_pixels(pixels2.begin(), pixels2.end());
  partitioner.finalize();

  std::size_t num_interactions = 0;
  for (const auto& [bid, _] : partitioner.block_index()) {
    const auto blk = partitioner.merge_blocks(bid);
    num_interactions += static_cast<std::size_t>(blk.nRecords);
  }

  CHECK(num_interactions == pixels1.size() + pixels2.size());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: SerializedBlockPQueue", "[hic][v9][short]") {
  using PQueue = SerializedBlockPQueue<std::uint64_t>;
  spdlog::set_level(spdlog::level::trace);
  const std::size_t num_threads = 32;

  std::vector<PQueue::Record> records(num_threads * 100);
  std::vector<std::uint64_t> blk_ids(num_threads * 100);

  std::random_device rd;
  std::mt19937_64 rand_eng(rd());

  std::uint64_t bid = 0;
  for (std::size_t i = 0; i < records.size(); ++i) {
    bid = std::uniform_int_distribution<std::uint64_t>{bid + 1, bid + 10}(rand_eng);
    records[i] = {bid, std::to_string(i), PQueue::Record::Status::SUCCESS};
    blk_ids[i] = bid;
  }

  PQueue queue(blk_ids.begin(), blk_ids.end(), num_threads - 1);

  REQUIRE(queue.size() == 0);
  REQUIRE(queue.capacity() > 0);

  std::atomic<std::size_t> i{};
  std::atomic<std::size_t> threads_started{};

  BS::thread_pool tpool(conditional_static_cast<BS::concurrency_t>(num_threads));

  auto producer = [&]() {
    std::random_device rd_;
    std::mt19937_64 rand_eng_(rd_());

    ++threads_started;
    while (threads_started != num_threads);

    while (true) {
      const auto idx = i++;
      if (idx >= records.size()) {
        return;
      }

      // Simulate time required for block compression
      const std::chrono::milliseconds sleep_time{
          std::uniform_int_distribution<std::int32_t>{25, 50}(rand_eng_)};
      std::this_thread::sleep_for(sleep_time);

      const auto& record = records[idx];
      while (!queue.try_enqueue(record.bid, record.serialized_block));
    }
  };

  std::vector<std::future<void>> producers(num_threads - 1);
  for (auto& prod : producers) {
    prod = tpool.submit_task(producer);
  }

  auto consumer = tpool.submit_task([&]() {
    std::vector<PQueue::Record> output;
    ++threads_started;
    while (true) {
      auto record = queue.dequeue_timed();
      if (record.status == PQueue::Record::Status::TIMEOUT) {
        continue;
      }
      if (record.status == PQueue::Record::Status::QUEUE_IS_CLOSED) {
        return output;
      }
      assert(!!record);
      output.emplace_back(std::move(record));
    }
  });

  for (auto& prod : producers) {
    prod.get();
  }
  const auto output = consumer.get();

  REQUIRE(output.size() == records.size());
  REQUIRE(i >= records.size());
  CHECK(queue.dequeue_timed().status == PQueue::Record::Status::QUEUE_IS_CLOSED);

  for (std::size_t j = 0; j < records.size(); ++j) {
    CHECK(output[j].bid == records[j].bid);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[maybe_unused]] static void compare_weights(const balancing::Weights& expected,
                                             const balancing::Weights& found) {
  REQUIRE(expected.size() == found.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (std::isnan(expected[i])) {
      CHECK(std::isnan(found[i]));
    } else {
      CHECK(expected[i] == found[i]);
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[maybe_unused]] static void hic_file_writer_compare_pixels(
    const std::vector<Pixel<float>>& expected, const std::vector<Pixel<float>>& found) {
  REQUIRE(expected.size() == found.size());

  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(expected[i].coords == found[i].coords);
    if (std::isnan(expected[i].count)) {
      CHECK(std::isnan(found[i].count));
    } else {
      CHECK_THAT(expected[i].count, Catch::Matchers::WithinRel(found[i].count));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[maybe_unused]] static void hic_file_writer_create_file_test(
    const std::string& path1, const std::string& path2,
    const std::vector<std::uint32_t>& resolutions, std::size_t num_threads,
    bool skip_all_vs_all_matrix) {
  {
    const auto chromosomes = hic::File(path1, resolutions.front()).chromosomes();
    const auto tmpdir = testdir() / (path1 + ".tmp");
    std::filesystem::create_directories(tmpdir);
    std::filesystem::remove(path2);
    HiCFileWriter w(path2, chromosomes, resolutions, "dm6", num_threads, 99'999, tmpdir, 1,
                    skip_all_vs_all_matrix);
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
      if (i % 2 == 0) {
        const auto resolution = resolutions[i];
        const hic::File f((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
        const auto sel1 = f.fetch("chr3R");
        const auto sel2 = f.fetch("chr3R", "chr4");
        w.add_pixels(resolution, sel1.begin<float>(), sel1.end<float>());
        w.add_pixels(resolution, sel2.begin<float>(), sel2.end<float>());
      }
    }
    w.serialize();
  }

  for (const auto& resolution : resolutions) {
    fmt::print(FMT_STRING("Comparing {}...\n"), resolution);
    const hic::File f1(path1, resolution);
    const hic::File f2(path2, resolution);

    const auto correct_pixels1 = f1.fetch("chr3R").read_all<float>();
    const auto correct_pixels2 = f1.fetch("chr3R", "chr4").read_all<float>();
    const auto pixels1 = f2.fetch("chr3R").read_all<float>();
    const auto pixels2 = f2.fetch("chr3R", "chr4").read_all<float>();

    hic_file_writer_compare_pixels(correct_pixels1, pixels1);
    hic_file_writer_compare_pixels(correct_pixels2, pixels2);

    const hic::File f3(path1, resolution, MatrixType::expected);
    const hic::File f4(path2, resolution, MatrixType::expected);

    const auto correct_expected_pixels1 = f3.fetch("chr3R").read_all<float>();
    const auto correct_expected_pixels2 = f4.fetch("chr3R", "chr4").read_all<float>();
    const auto expected_pixels1 = f3.fetch("chr3R").read_all<float>();
    const auto expected_pixels2 = f4.fetch("chr3R", "chr4").read_all<float>();

    // NOLINTNEXTLINE(*-suspicious-call-argument)
    hic_file_writer_compare_pixels(correct_expected_pixels1, expected_pixels1);
    // NOLINTNEXTLINE(*-suspicious-call-argument)
    hic_file_writer_compare_pixels(correct_expected_pixels2, expected_pixels2);
  }
}

[[nodiscard]] static std::vector<double> generate_random_weights(std::uint32_t chrom_size,
                                                                 std::uint32_t resolution) {
  std::vector<double> buff((chrom_size + resolution - 1) / resolution);
  std::random_device rd{};
  std::mt19937_64 rand_eng{rd()};

  std::generate(buff.begin(), buff.end(),
                [&]() { return std::uniform_real_distribution<double>{0, 1}(rand_eng); });

  return buff;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: HiCFileWriter (creation)", "[hic][v9][long]") {
  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_writer_001.hic").string();
  const auto path3 = (testdir() / "hic_writer_002.hic").string();

  spdlog::set_level(spdlog::level::trace);

  SECTION("create file (mt)") {
    const std::vector<std::uint32_t> resolutions{25'000, 1'000'000, 2'500'000};
    hic_file_writer_create_file_test(path1, path2, resolutions, 16, true);
  }

  SECTION("regression PR 180") {
    // Ensure we can create .hic files having bin tables with 1 bin per chromosome
    // See https://github.com/paulsengroup/hictk/pull/180
    const hictk::Reference chromosomes{{0, "chr1", 10}};
    HiCFileWriter w(path3, chromosomes, {100});

    const std::vector<Pixel<float>> pixels{Pixel<float>{w.bins(100), 0, 0, 1.0F}};
    w.add_pixels(100, pixels.begin(), pixels.end());
    w.serialize();  // Before PR 180, this used to throw
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: HiCFileWriter (add weights)", "[hic][v9][long]") {
  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_writer_004.hic").string();

  SECTION("add weights") {
    const std::uint32_t resolution = 500'000;
    const hic::File hf1(path1, resolution);

    {
      // init file
      HiCFileWriter w(path2, hf1.chromosomes(), {hf1.resolution()}, "dm6");
      const auto sel = hf1.fetch();
      w.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
      w.serialize();
    }

    // add normalization weights
    {
      HiCFileWriter w(path2);
      for (const auto& chrom : w.chromosomes()) {
        if (chrom.is_all()) {
          continue;
        }
        w.add_norm_vector("SCALE", chrom, "BP", hf1.resolution(),
                          hf1.normalization("SCALE", chrom));
      }

      CHECK_THROWS_WITH(
          w.add_norm_vector("SCALE", hf1.chromosomes().at("chr2L"), "BP", hf1.resolution(),
                            hf1.normalization("SCALE", hf1.chromosomes().at("chr2L"))),
          Catch::Matchers::ContainsSubstring("file already contains"));

      CHECK_THROWS_WITH(
          w.add_norm_vector("VC", w.chromosomes().at("chr2L"), "BP", hf1.resolution(),
                            balancing::Weights{{1, 2, 3}, balancing::Weights::Type::DIVISIVE}),
          Catch::Matchers::ContainsSubstring("weight shape mismatch"));

      w.write_norm_vectors_and_norm_expected_values();
    }

    // compare
    const hic::File hf2(path2, resolution);

    const auto avail_norms = hf2.avail_normalizations();
    REQUIRE(avail_norms.size() == 1);
    CHECK(avail_norms.front() == balancing::Method::SCALE());

    compare_weights(hf1.normalization("SCALE"), hf2.normalization("SCALE"));
    hic_file_writer_compare_pixels(hf1.fetch(balancing::Method::SCALE()).read_all<float>(),
                                   hf2.fetch(balancing::Method::SCALE()).read_all<float>());

    const hic::File hf3(path1, resolution, MatrixType::expected);
    const hic::File hf4(path2, resolution, MatrixType::expected);

    compare_weights(hf3.normalization("SCALE"), hf4.normalization("SCALE"));
    hic_file_writer_compare_pixels(hf3.fetch(balancing::Method::SCALE()).read_all<float>(),
                                   hf4.fetch(balancing::Method::SCALE()).read_all<float>());
  }

  SECTION("overwrite weights") {
    const std::uint32_t resolution = 500'000;

    {
      // init file
      const hic::File hf(path1, resolution);
      std::filesystem::remove(path2);  // NOLINT
      HiCFileWriter w(path2, hf.chromosomes(), {hf.resolution()}, "dm6");
      const auto sel = hf.fetch();
      w.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
      w.serialize();
    }

    // add normalization weights
    std::vector<double> weights{};
    {
      HiCFileWriter w(path2);
      for (const auto& chrom : w.chromosomes()) {
        if (chrom.is_all()) {
          continue;
        }
        const auto buff = generate_random_weights(chrom.size(), resolution);
        weights.insert(weights.end(), buff.begin(), buff.end());
        w.add_norm_vector("FOO", chrom, "BP", resolution,
                          balancing::Weights{buff, balancing::Weights::Type::DIVISIVE});
      }
    }
    // compare weights
    {
      const hic::File hf(path2, resolution);
      compare_weights(balancing::Weights{weights, balancing::Weights::Type::DIVISIVE},
                      hf.normalization("FOO"));
    }

    // overwrite weights
    {
      weights.clear();
      HiCFileWriter w(path2);
      for (const auto& chrom : w.chromosomes()) {
        if (chrom.is_all()) {
          continue;
        }
        const auto buff = generate_random_weights(chrom.size(), resolution);
        weights.insert(weights.end(), buff.begin(), buff.end());
        w.add_norm_vector("FOO", chrom, "BP", resolution,
                          balancing::Weights{buff, balancing::Weights::Type::DIVISIVE});
      }
    }

    // compare weights
    const hic::File hf(path2, resolution);
    compare_weights(balancing::Weights{weights, balancing::Weights::Type::DIVISIVE},
                    hf.normalization("FOO"));
  }
}

}  // namespace hictk::hic::test::file_writer
