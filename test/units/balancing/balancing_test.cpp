// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
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
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "hictk/balancing/ice.hpp"
#include "hictk/balancing/scale.hpp"
#include "hictk/balancing/sparse_matrix.hpp"
#include "hictk/balancing/vc.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/file.hpp"
#include "hictk/filestream.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "tmpdir.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::balancing {

[[nodiscard]] static hictk::balancing::Weights read_weights(const std::filesystem::path& path,
                                                            hictk::balancing::Weights::Type type,
                                                            char sep = '\n') {
  assert(std::filesystem::exists(path));
  std::ifstream ifs(path);
  std::string strbuf;
  std::vector<double> buffer{};

  while (std::getline(ifs, strbuf, sep)) {
    buffer.push_back(std::stod(strbuf));
  }

  return {buffer, type};
}

static void compare_weights(const hictk::balancing::Weights& weights_,
                            const hictk::balancing::Weights& expected_, double tol = 5.0e-3) {
  REQUIRE(weights_.size() == expected_.size());

  const auto weights = weights_(hictk::balancing::Weights::Type::DIVISIVE);
  const auto expected = expected_(hictk::balancing::Weights::Type::DIVISIVE);

  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (std::isnan(expected[i])) {
      CHECK(std::isnan(weights[i]));
    } else {
      CHECK_THAT(weights[i], Catch::Matchers::WithinRel(expected[i], tol));
    }
  }
}

template <typename T>
static void compare_vectors(const std::vector<T>& v1, const std::vector<T>& v2) {
  REQUIRE(v1.size() == v2.size());

  for (std::size_t i = 0; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: VectorOfAtomicDecimals", "[balancing][short]") {
  using VectorOfAtomicDecimals = hictk::balancing::internal::VectorOfAtomicDecimals;

  SECTION("Ctors") {
    VectorOfAtomicDecimals v1(10);
    CHECK(v1.size() == 10);
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
    CHECK(v1.decimal_bits() == 22);

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
  }

  SECTION("atomic modifiers") {
    SECTION("add (st)") {
      VectorOfAtomicDecimals v(10);

      v.add(0, 0);
      CHECK(v[0] == 0);

      v.add(0, 1.0e-3);
      CHECK_THAT(v[0], Catch::Matchers::WithinAbs(1.0e-3, 1.0e-6));

      v.add(0, 100.0e9);
      CHECK_THAT(v[0], Catch::Matchers::WithinRel(100.0e9));

      v.add(0, v.domain(false).second - 100.0e9 + 1);
      CHECK(std::isinf(v[0]));

      v.add(0, std::numeric_limits<double>::quiet_NaN());
      CHECK(std::isnan(v[0]));

      v.add(0, 10);
      CHECK(std::isnan(v[0]));

      v.add(0, std::numeric_limits<double>::infinity());
      CHECK(std::isnan(v[0]));
    }

    SECTION("add (mt, wo/ overflow)") {
      VectorOfAtomicDecimals v(1);

      auto worker = [&](std::size_t num_threads, std::atomic<std::size_t>& threads_started,
                        std::size_t iters = 1'000'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        double tot = 0.0;

        ++threads_started;
        while (threads_started != num_threads);
        for (std::size_t i = 0; i < iters; ++i) {
          const auto n = std::uniform_real_distribution<double>{0, 10}(rand_gen);
          v.add(0, n);
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
                        std::size_t iters = 100'000) {
        std::random_device rd{};
        std::mt19937_64 rand_gen{rd()};

        const auto ub = v.domain(false).second / static_cast<double>(iters / num_threads);

        double tot = 0.0;

        ++threads_started;
        while (threads_started != num_threads);

        for (std::size_t i = 0; i < iters; ++i) {
          const auto n = std::uniform_real_distribution<double>{0, ub}(rand_gen);
          v.add(0, n);
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
      auto f = filestream::FileStream::create(tmpfile.string());

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

      std::filesystem::remove(tmpfile);
      auto f = filestream::FileStream::create(tmpfile.string());

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
  const auto tmpfile = testdir() / "sparse_matrix_chunked.tmp";

  SECTION("accessors") { CHECK(SparseMatrixChunked{tmpfile, 2, 0}.empty()); }

  SECTION("push_back") {
    SparseMatrixChunked m{tmpfile, 2, 0};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();

    CHECK(m.size() == pixels.size());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: ICE (intra)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_ice_intra.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.ICE.cis.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::ICE::Type::cis;
        const auto weights = hictk::balancing::ICE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }

      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::ICE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::ICE::Type::cis;
        const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::ICE(var_bin_file, hictk::balancing::ICE::Type::cis));
    CHECK_THROWS(hictk::balancing::ICE(storage_mode_square_file, hictk::balancing::ICE::Type::cis));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: ICE (inter)", "[balancing][medium]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_ice_inter.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.ICE.trans.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::ICE::Type::trans;
        const auto weights = hictk::balancing::ICE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }

      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::ICE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::ICE::Type::trans;
        const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::ICE(var_bin_file, hictk::balancing::ICE::Type::trans));
    CHECK_THROWS(
        hictk::balancing::ICE(storage_mode_square_file, hictk::balancing::ICE::Type::trans));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: ICE (gw)", "[balancing][medium]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_ice_gw.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.ICE.gw.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::ICE::Type::gw;
        const auto weights = hictk::balancing::ICE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }

      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::ICE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::ICE::Type::gw;
        const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::MULTIPLICATIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::ICE(var_bin_file, hictk::balancing::ICE::Type::gw));
    CHECK_THROWS(hictk::balancing::ICE(storage_mode_square_file, hictk::balancing::ICE::Type::gw));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: VC (intra)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.VC.cis.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      constexpr auto type = hictk::balancing::VC::Type::cis;
      const auto weights = hictk::balancing::VC(f, type).get_weights();
      const auto expected_weights =
          read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

      compare_weights(weights, expected_weights);
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::VC(var_bin_file, hictk::balancing::VC::Type::cis));
    CHECK_THROWS(hictk::balancing::VC(storage_mode_square_file, hictk::balancing::VC::Type::cis));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: VC (inter)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.VC.inter.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      constexpr auto type = hictk::balancing::VC::Type::trans;
      const auto weights = hictk::balancing::VC(f, type).get_weights();
      const auto expected_weights =
          read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

      compare_weights(weights, expected_weights);
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::VC(var_bin_file, hictk::balancing::VC::Type::trans));
    CHECK_THROWS(hictk::balancing::VC(storage_mode_square_file, hictk::balancing::VC::Type::trans));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: VC (gw)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.VC.gw.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      constexpr auto type = hictk::balancing::VC::Type::gw;
      const auto weights = hictk::balancing::VC(f, type).get_weights();
      const auto expected_weights =
          read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

      compare_weights(weights, expected_weights);
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::VC(var_bin_file, hictk::balancing::VC::Type::gw));
    CHECK_THROWS(hictk::balancing::VC(storage_mode_square_file, hictk::balancing::VC::Type::gw));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: SCALE (intra)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_scale_cis.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.SCALE.cis.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::SCALE::Type::cis;
        const auto weights = hictk::balancing::SCALE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }
      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::SCALE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::SCALE::Type::cis;
        const auto weights = hictk::balancing::SCALE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::SCALE(var_bin_file, hictk::balancing::SCALE::Type::cis));
    CHECK_THROWS(
        hictk::balancing::SCALE(storage_mode_square_file, hictk::balancing::SCALE::Type::cis));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: SCALE (inter)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_scale_trans.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.SCALE.inter.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::SCALE::Type::trans;
        const auto weights = hictk::balancing::SCALE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }

      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::SCALE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::SCALE::Type::trans;
        const auto weights = hictk::balancing::SCALE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::SCALE(var_bin_file, hictk::balancing::SCALE::Type::trans));
    CHECK_THROWS(
        hictk::balancing::SCALE(storage_mode_square_file, hictk::balancing::SCALE::Type::trans));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: SCALE (gw)", "[balancing][short]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.2500000.hic")};

  const auto tmpfile = testdir() / "balancing_scale_gw.tmp";
  const auto path_weights = datadir / "balancing/ENCFF993FGR.2500000.SCALE.gw.txt";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        constexpr auto type = hictk::balancing::SCALE::Type::gw;
        const auto weights = hictk::balancing::SCALE(f, type).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }

      SECTION("chunked") {
        std::filesystem::remove(tmpfile);
        auto params = hictk::balancing::SCALE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        constexpr auto type = hictk::balancing::SCALE::Type::gw;
        const auto weights = hictk::balancing::SCALE(f, type, params).get_weights();
        const auto expected_weights =
            read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

        compare_weights(weights, expected_weights);
      }
    }
  }

  SECTION("invalid files") {
    const cooler::File var_bin_file(
        (datadir / "cooler/cooler_variable_bins_test_file.cool").string());
    const cooler::File storage_mode_square_file(
        (datadir / "cooler/cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());

    CHECK_THROWS(hictk::balancing::SCALE(var_bin_file, hictk::balancing::SCALE::Type::gw));
    CHECK_THROWS(
        hictk::balancing::SCALE(storage_mode_square_file, hictk::balancing::SCALE::Type::gw));
  }
}

TEST_CASE("Balancing: SCALE (edge cases)", "[balancing][medium]") {
  SECTION("diverged") {
    const auto path = datadir / "hic/4DNFIZ1ZVXC8.hic9";
    const auto path_weights = datadir / "balancing/4DNFIZ1ZVXC8.chr2L.10000.SCALE.txt";

    const hictk::File f(path.string(), 10'000);
    const auto sel = f.fetch("chr2L");
    const auto weights =
        hictk::balancing::SCALE(sel.template begin<double>(), sel.template end<double>(),
                                f.bins().subset("chr2L"))
            .get_weights();
    const auto expected_weights =
        read_weights(path_weights, hictk::balancing::Weights::Type::DIVISIVE);

    compare_weights(weights, expected_weights);
  }
}

}  // namespace hictk::test::balancing
