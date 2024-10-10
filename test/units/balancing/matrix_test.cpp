// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <zstd.h>

#include <array>
#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <future>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "./common.hpp"
#include "hictk/balancing/sparse_matrix.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/filestream.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "tmpdir.hpp"

namespace hictk::test::balancing {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Balancing: AtomicBitSet", "[balancing][short]") {
  using AtomicBitSet = hictk::balancing::internal::AtomicBitSet;

  SECTION("Ctors") {
    AtomicBitSet b1(10);
    REQUIRE(b1.size() == 10);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK_FALSE(b1.atomic_test(i));
    }

    b1.atomic_set(0, true);
    CHECK(b1.atomic_test(0));

    const AtomicBitSet b2(b1);
    REQUIRE(b2.size() == b1.size());
    CHECK(b2.atomic_test(0));
    CHECK_FALSE(b2.atomic_test(1));
  }

  SECTION("operator=") {
    AtomicBitSet b1(10);
    b1.atomic_set(0, true);

    REQUIRE(b1.atomic_test(0));

    const AtomicBitSet b2(1);
    REQUIRE_FALSE(b2.atomic_test(0));

    b1 = b2;
    CHECK(b1.size() == b2.size());
    CHECK_FALSE(b1.atomic_test(0));
  }

  SECTION("accessors") {
    const AtomicBitSet b(10, true);

    CHECK(b.size() == 10);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(b.atomic_test(i));
    }
  }

  SECTION("non-atomic modifiers") {
    AtomicBitSet b(10);
    b.resize(15, true);

    for (std::size_t i = 0; i < 15; ++i) {
      if (i < 10) {
        CHECK_FALSE(b.atomic_test(i));
      } else {
        CHECK(b.atomic_test(i));
      }
    }

    b.fill(false);
    for (std::size_t i = 0; i < 15; ++i) {
      CHECK_FALSE(b.atomic_test(i));
    }
  }

  SECTION("atomic modifiers") {
    constexpr std::size_t nthreads{2};

    SECTION("concurrent set") {
      std::atomic<bool> buff{false};
      AtomicBitSet b(1);

      auto worker = [&](std::size_t num_threads, std::atomic<std::size_t>& threads_started,
                        std::size_t iters = 1'000'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        ++threads_started;
        // clang-format off
        while (threads_started != num_threads);  // NOLINT
        // clang-format on
        for (std::size_t i = 0; i < iters; ++i) {
          const auto x = std::bernoulli_distribution{}(rand_gen);
          b.atomic_set(0, x);
          buff = x;
        }
      };

      for (std::size_t i = 0; i < 10; ++i) {
        b.fill(false);
        buff = false;

        std::atomic<std::size_t> threads_started{0};
        std::vector<std::future<void>> futures(nthreads);
        std::generate(futures.begin(), futures.end(),
                      [&]() { return std::async(worker, nthreads, std::ref(threads_started)); });

        std::for_each(futures.begin(), futures.end(), [](auto& fut) { fut.get(); });

        CHECK(buff.load() == b.atomic_test(0));
      }
    }

    SECTION("concurrent set adjacent") {
      std::vector<std::atomic<bool>> buff(4);
      AtomicBitSet b(4);

      auto worker = [&](std::size_t num_threads, std::atomic<std::size_t>& threads_started,
                        std::size_t iters = 1'000'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        ++threads_started;
        // clang-format off
        while (threads_started != num_threads);  // NOLINT
        // clang-format on
        for (std::size_t i = 0; i < iters; ++i) {
          const auto j = std::uniform_int_distribution<std::size_t>{0, b.size() - 1}(rand_gen);
          const auto x = std::bernoulli_distribution{}(rand_gen);
          b.atomic_set(j, x);
          buff[j] = x;
        }
      };

      for (std::size_t i = 0; i < 10; ++i) {
        b.fill(false);
        std::fill(buff.begin(), buff.end(), false);

        std::atomic<std::size_t> threads_started{0};
        std::vector<std::future<void>> futures(nthreads);
        std::generate(futures.begin(), futures.end(),
                      [&]() { return std::async(worker, nthreads, std::ref(threads_started)); });

        std::for_each(futures.begin(), futures.end(), [](auto& fut) { fut.get(); });

        for (std::size_t j = 0; j < b.size(); ++j) {
          CHECK(buff[j].load() == b.atomic_test(j));
        }
      }
    }
  }
}

TEST_CASE("Balancing: VectorOfAtomicDecimals", "[balancing][short]") {
  using VectorOfAtomicDecimals = hictk::balancing::internal::VectorOfAtomicDecimals;

  SECTION("Ctors") {
    VectorOfAtomicDecimals v1(10);
    REQUIRE(v1.size() == 10);
    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == 0);
    }

    v1.set(0, 10);
    REQUIRE_THAT(v1[0], Catch::Matchers::WithinRel(10.0));

    const VectorOfAtomicDecimals v2(v1);
    CHECK(v2.size() == 10);
    REQUIRE_THAT(v2[0], Catch::Matchers::WithinRel(10.0));
  }

  SECTION("operator=") {
    VectorOfAtomicDecimals v1(10);
    v1.set(0, 10);
    REQUIRE_THAT(v1[0], Catch::Matchers::WithinRel(10.0));

    VectorOfAtomicDecimals v2(1);
    REQUIRE(v2.size() == 1);
    v2 = v1;

    CHECK(v2.size() == v1.size());
    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i]);
    }
  }

  SECTION("accessors") {
    VectorOfAtomicDecimals v1(10);
    v1.set(0, 10);
    REQUIRE_THAT(v1[0], Catch::Matchers::WithinRel(10.0));

    CHECK(v1.size() == 10);
    CHECK(!v1.empty());
    CHECK(v1.decimal_bits() == 30);

    const auto& v2 = v1();
    REQUIRE(v1.size() == v2.size());
    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK_THAT(v1[i], Catch::Matchers::WithinRel(v2[i]));
    }
  }

  SECTION("non-atomic modifiers") {
    VectorOfAtomicDecimals v1(10);
    for (std::size_t i = 0; i < v1.size(); ++i) {
      v1.set(i, static_cast<double>(i));
    }

    SECTION("resize") {
      v1.resize(20);
      REQUIRE(v1.size() == 20);

      for (std::size_t i = 0; i < v1.size(); ++i) {
        if (i < 10) {
          CHECK_THAT(v1[i], Catch::Matchers::WithinRel(static_cast<double>(i)));
        } else {
          CHECK(v1[i] == 0);
        }
      }

      v1.resize(5);
      REQUIRE(v1.size() == 5);
      for (std::size_t i = 0; i < v1.size(); ++i) {
        CHECK_THAT(v1[i], Catch::Matchers::WithinRel(static_cast<double>(i)));
      }
    }

    SECTION("fill") {
      v1.fill(17);
      REQUIRE(v1.size() == 10);
      for (std::size_t i = 0; i < v1.size(); ++i) {
        CHECK_THAT(v1[i], Catch::Matchers::WithinRel(17.0));
      }
    }

    SECTION("multiply") {
      SECTION("finite") {
        const std::vector<double> vfinite1(v1.size(), 10);
        v1.fill(17);
        v1.multiply(vfinite1);

        REQUIRE(v1.size() == 10);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK_THAT(v1[i], Catch::Matchers::WithinRel(170.0));
        }

        const std::vector<double> vfinite2(v1.size(), 0);
        v1.fill(17);
        v1.multiply(vfinite2);

        REQUIRE(v1.size() == 10);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK(v1[i] == 0.0);
        }

        const auto max_value = v1.domain(false).second;
        const std::vector<double> vfinite3(v1.size(), max_value);

        v1.fill(1);
        v1.multiply(vfinite3);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK_THAT(v1[i], Catch::Matchers::WithinRel(max_value));
        }

        const std::vector<double> vfinite4(v1.size(), std::nextafter(max_value, max_value + 1));

        v1.fill(1);
        v1.multiply(vfinite4);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK(std::isinf(v1[i]));
        }
      }

      SECTION("nan") {
        const std::vector<double> vnan(v1.size(), std::numeric_limits<double>::quiet_NaN());

        v1.fill(17);
        v1.multiply(vnan);
        REQUIRE(v1.size() == 10);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK(std::isnan(v1[i]));
        }
      }

      SECTION("inf") {
        const std::vector<double> vinf(v1.size(), std::numeric_limits<double>::infinity());

        v1.fill(17);
        v1.multiply(vinf);
        REQUIRE(v1.size() == 10);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK(std::isinf(v1[i]));
        }

        v1.fill(0);
        v1.multiply(vinf);
        REQUIRE(v1.size() == 10);
        for (std::size_t i = 0; i < v1.size(); ++i) {
          CHECK(std::isnan(v1[i]));
        }
      }
    }

    SECTION("set") {
      VectorOfAtomicDecimals v(10);

      v.set(0, 0);
      CHECK(v[0] == 0);

      v.set(0, 1.0e-3);
      CHECK_THAT(v[0], Catch::Matchers::WithinAbs(1.0e-3, 1.0e-6));

      v.set(0, 1.0e9);
      CHECK_THAT(v[0], Catch::Matchers::WithinRel(1.0e9));

      v.set(0, v.domain(false).second + 1);
      CHECK(std::isinf(v[0]));

      v.set(0, std::numeric_limits<double>::quiet_NaN());
      CHECK(std::isnan(v[0]));

      v.set(0, std::numeric_limits<double>::infinity());
      CHECK(std::isinf(v[0]));

      v.set(0, 0);
      CHECK(v[0] == 0);
    }
  }

  SECTION("atomic modifiers") {
    SECTION("add (st)") {
      VectorOfAtomicDecimals v(10);

      v.atomic_add(0, 0);
      CHECK(v[0] == 0);

      v.atomic_add(0, 1.0e-3);
      CHECK_THAT(v[0], Catch::Matchers::WithinAbs(1.0e-3, 1.0e-6));

      v.set(0, 0);
      v.atomic_add(0, 10.0e9);
      CHECK_THAT(v[0], Catch::Matchers::WithinAbs(10.0e9, 1.0e-6));

      v.atomic_add(0, v.domain(false).second - 10.0e9 + 1000);
      CHECK(std::isinf(v[0]));

      v.atomic_add(0, std::numeric_limits<double>::quiet_NaN());
      CHECK(std::isnan(v[0]));

      v.atomic_add(0, 10);
      CHECK(std::isnan(v[0]));

      v.atomic_add(0, std::numeric_limits<double>::infinity());
      CHECK(std::isnan(v[0]));
    }

    SECTION("add (mt, wo/ overflow)") {
      VectorOfAtomicDecimals v(1);

      auto worker = [&](std::size_t num_threads, std::atomic<std::size_t>& threads_started,
                        std::size_t iters = 100'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        double tot = 0.0;

        ++threads_started;
        // clang-format off
        while (threads_started != num_threads);  // NOLINT
        // clang-format on
        for (std::size_t i = 0; i < iters; ++i) {
          const auto n = std::uniform_real_distribution<double>{0, 10}(rand_gen);
          v.atomic_add(0, n);
          tot += n;
        }

        return tot;
      };

      for (std::size_t i = 0; i < 10; ++i) {
        v.fill(0);

        std::vector<std::future<double>> futures(2);
        std::atomic<std::size_t> threads_started{};

        std::generate(futures.begin(), futures.end(), [&]() {
          return std::async(worker, futures.size(), std::ref(threads_started));
        });
        const auto tot =
            std::accumulate(futures.begin(), futures.end(), 0.0,
                            [&](double accumulator, auto& fut) { return accumulator + fut.get(); });

        REQUIRE(tot <= v.domain(false).second);
        CHECK_THAT(v[0], Catch::Matchers::WithinAbs(tot, 1.0e5));
      }
    }

    SECTION("add (mt, w/ overflow)") {
      VectorOfAtomicDecimals v(1);

      auto worker = [&](std::size_t num_threads, std::atomic<std::size_t>& threads_started,
                        std::size_t iters = 1'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        // NOLINTNEXTLINE(bugprone-integer-division)
        const auto ub = v.domain(false).second / static_cast<double>(iters / num_threads);

        double tot = 0.0;

        ++threads_started;
        // clang-format off
        while (threads_started != num_threads);  // NOLINT
        // clang-format on

        for (std::size_t i = 0; i < iters; ++i) {
          const auto n = std::uniform_real_distribution<double>{0, ub}(rand_gen);
          v.atomic_add(0, n);
          tot += n;
        }

        return tot;
      };

      for (std::size_t i = 0; i < 100; ++i) {
        v.fill(0);

        std::vector<std::future<double>> futures(2);
        std::atomic<std::size_t> threads_started{};

        std::generate(futures.begin(), futures.end(), [&]() {
          return std::async(worker, futures.size(), std::ref(threads_started));
        });
        const auto tot =
            std::accumulate(futures.begin(), futures.end(), 0.0,
                            [&](double accumulator, auto& fut) { return accumulator + fut.get(); });

        if (tot > v.domain(false).second) {
          CHECK(std::isinf(v[0]));
        } else {
          CHECK_THAT(v[0], Catch::Matchers::WithinAbs(tot, 1.0e-5));
        }
      }
    }
  }

  SECTION("encode/decode") {
    std::random_device rd{};
    std::mt19937_64 rand_eng{rd()};

    VectorOfAtomicDecimals v(1);

    SECTION("small numbers") {
      std::uniform_real_distribution<double> dist{0.0, 10.0};
      for (std::size_t i = 0; i < 500'000; ++i) {
        const auto n = dist(rand_eng);
        v.set(0, n);
        CHECK_THAT(v[0], Catch::Matchers::WithinAbs(n, 1.0e6));
      }
    }

    SECTION("intermediate numbers") {
      std::uniform_real_distribution<double> dist{10.0, 1.0e6};
      for (std::size_t i = 0; i < 500'000; ++i) {
        const auto n = dist(rand_eng);
        v.set(0, n);
        CHECK_THAT(v[0], Catch::Matchers::WithinAbs(n, 1.0e6));
      }
    }

    SECTION("large numbers") {
      const auto ub = v.domain(false).second;
      // Draw a number outside the represented domain ~10% of the times
      std::uniform_real_distribution<double> dist{1.0e6, ub * 1.1};

      for (std::size_t i = 0; i < 500'000; ++i) {
        const auto n = dist(rand_eng);
        v.set(0, n);
        if (n > ub) {
          CHECK(std::isinf(v[0]));
        } else {
          CHECK_THAT(v[0], Catch::Matchers::WithinAbs(n, 1.0e6));
        }
      }
    }
  }
}

TEST_CASE("Balancing: SparseMatrix", "[balancing][short]") {
  using SparseMatrix = hictk::balancing::internal::SparseMatrix;
  const BinTable bins{Reference{Chromosome{0, "chr0", 50}, Chromosome{1, "chr1", 100},
                                Chromosome{2, "chr2", 50}, Chromosome{3, "chr3", 50}},
                      50};
  // clang-format off
const std::vector<ThinPixel<std::int32_t>> pixels{
{1, 1, 1}, {1, 2, 2}, {2, 2, 3},  // chr1
{3, 3, 4}, {3, 4, 5}};            // chr2
  // clang-format on

  SECTION("accessors") { CHECK(SparseMatrix{}.empty()); }

  SECTION("push_back") {
    SparseMatrix m{};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();
    CHECK(m.size() == pixels.size());

    m.clear();
    CHECK(m.empty());
  }

  SECTION("serde") {
    const auto tmpfile = testdir() / "sparse_matrix_serde.bin";
    std::unique_ptr<ZSTD_CCtx_s> zstd_cctx{ZSTD_createCCtx()};
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx{ZSTD_createDCtx()};

    std::string buff{};

    SECTION("empty matrix") {
      auto f = filestream::FileStream<>::create(tmpfile.string(), nullptr);

      SparseMatrix m1{};
      SparseMatrix m2{};
      m1.finalize();
      m1.serialize(f, buff, *zstd_cctx);
      f.seekg(std::ios::beg);
      m2.deserialize(f, buff, *zstd_dctx);

      compare_vectors(m1.bin1_ids(), m2.bin1_ids());
      compare_vectors(m1.bin2_ids(), m2.bin2_ids());
      compare_vectors(m1.counts(), m2.counts());
    }

    SECTION("full matrix") {
      SparseMatrix m1{};
      for (const auto& p : pixels) {
        m1.push_back(p.bin1_id, p.bin2_id, p.count);
      }
      m1.finalize();

      std::filesystem::remove(tmpfile);  // NOLINT
      auto f = filestream::FileStream<>::create(tmpfile.string(), nullptr);

      SparseMatrix m2{};
      m1.serialize(f, buff, *zstd_cctx);
      f.seekg(std::ios::beg);
      m2.deserialize(f, buff, *zstd_dctx);

      compare_vectors(m1.bin1_ids(), m2.bin1_ids());
      compare_vectors(m1.bin2_ids(), m2.bin2_ids());
      compare_vectors(m1.counts(), m2.counts());
    }
  }
}

TEST_CASE("Balancing: SparseMatrixChunked", "[balancing][short]") {
  using SparseMatrixChunked = hictk::balancing::internal::SparseMatrixChunked;
  const BinTable bins{Reference{Chromosome{0, "chr0", 50}, Chromosome{1, "chr1", 100},
                                Chromosome{2, "chr2", 50}, Chromosome{3, "chr3", 50}},
                      50};
  // clang-format off
const std::vector<ThinPixel<std::int32_t>> pixels{
{1, 1, 1}, {1, 2, 2}, {2, 2, 3},  // chr1
{3, 3, 4}, {3, 4, 5}};            // chr2
  // clang-format on

  SECTION("accessors") {
    CHECK(SparseMatrixChunked{}.empty());
    CHECK(SparseMatrixChunked{}.num_chunks() == 0);
  }

  SECTION("push_back") {
    SparseMatrixChunked m{2};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();
    CHECK(m.size() == pixels.size());
    CHECK(m.num_chunks() == (m.size() + m.chunk_size() - 1) / m.chunk_size());

    m.clear();
    CHECK(m.empty());
    CHECK(m.num_chunks() == 0);
  }
}

TEST_CASE("Balancing: FileBackedSparseMatrix", "[balancing][short]") {
  using FileBackedSparseMatrix = hictk::balancing::internal::FileBackedSparseMatrix;
  const BinTable bins{Reference{Chromosome{0, "chr0", 50}, Chromosome{1, "chr1", 100},
                                Chromosome{2, "chr2", 50}, Chromosome{3, "chr3", 50}},
                      50};
  // clang-format off
const std::vector<ThinPixel<std::int32_t>> pixels{
{1, 1, 1}, {1, 2, 2}, {2, 2, 3},  // chr1
{3, 3, 4}, {3, 4, 5}};            // chr2
  // clang-format on
  const auto tmpfile = testdir() / "sparse_matrix_chunked.tmp";

  SECTION("accessors") { CHECK(FileBackedSparseMatrix{tmpfile, 2, 0}.empty()); }

  SECTION("push_back") {
    FileBackedSparseMatrix m{tmpfile, 2, 0};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();

    CHECK(m.size() == pixels.size());
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::balancing
