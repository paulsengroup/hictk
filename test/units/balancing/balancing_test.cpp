// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>

#include "hictk/balancing/ice.hpp"
#include "hictk/balancing/methods.hpp"
#include "hictk/file.hpp"
#include "tmpdir.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::balancing {

[[nodiscard]] static std::vector<double> read_weights(const std::filesystem::path& path,
                                                      char sep = '\n') {
  assert(std::filesystem::exists(path));
  std::ifstream ifs(path);
  std::string strbuf;
  std::vector<double> buffer{};

  while (std::getline(ifs, strbuf, sep)) {
    buffer.push_back(std::stod(strbuf));
  }

  return buffer;
}

static void compare_weights(const std::vector<double>& weights, const std::vector<double>& expected,
                            double tol = 1.0e-6) {
  REQUIRE(weights.size() == expected.size());

  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (std::isnan(weights[i])) {
      CHECK(std::isnan(expected[i]));
    } else {
      CHECK_THAT(weights[i], Catch::Matchers::WithinAbs(expected[i], tol));
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
TEST_CASE("Balancing: SparseMatrix") {
  using SparseMatrix = hictk::balancing::SparseMatrix;
  const BinTable bins{Reference{Chromosome{0, "chr0", 50}, Chromosome{1, "chr1", 100},
                                Chromosome{2, "chr2", 50}, Chromosome{3, "chr3", 50}},
                      50};
  // clang-format off
  const std::vector<ThinPixel<std::int32_t>> pixels{
      {1, 0, 1}, {1, 1, 2}, {2, 1, 3},  // chr1
      {3, 0, 4}, {3, 1, 5}};            // chr2
  // clang-format on

  SECTION("accessors") {
    CHECK(SparseMatrix{}.empty());
    CHECK(SparseMatrix{bins}.empty());
  }

  SECTION("push_back") {
    SparseMatrix m{bins};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();
    CHECK(m.size() == pixels.size());

    m.clear();
    CHECK(m.empty());
  }

  SECTION("subset") {
    SparseMatrix m{bins};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();

    CHECK(m.subset(0).empty());
    CHECK(m.subset(1).size() == 3);
    CHECK(m.subset(2).size() == 2);
    CHECK(m.subset(3).empty());
  }

  SECTION("serde") {
    const auto tmpfile = testdir() / "sparse_matrix_serde.bin";

    SECTION("empty matrix") {
      std::fstream f{};
      f.open(tmpfile, std::ios::in | std::ios::out | std::ios::trunc);
      f.exceptions(std::ios::badbit | std::ios::failbit);

      SparseMatrix m1{};
      SparseMatrix m2{};
      m1.finalize();
      m1.serialize(f);
      f.seekg(std::ios::beg);
      m2.deserialize(f);

      compare_vectors(m1.bin1_ids(), m2.bin1_ids());
      compare_vectors(m1.bin2_ids(), m2.bin2_ids());
      compare_vectors(m1.counts(), m2.counts());
      compare_vectors(m1.chrom_offsets(), m2.chrom_offsets());
    }

    SECTION("full matrix") {
      SparseMatrix m1{bins};
      for (const auto& p : pixels) {
        m1.push_back(p.bin1_id, p.bin2_id, p.count);
      }
      m1.finalize();

      std::fstream f{};
      f.open(tmpfile, std::ios::in | std::ios::out | std::ios::trunc);
      f.exceptions(std::ios::badbit | std::ios::failbit);

      SparseMatrix m2{bins};
      m1.serialize(f);
      f.seekg(std::ios::beg);
      m2.deserialize(f);

      compare_vectors(m1.bin1_ids(), m2.bin1_ids());
      compare_vectors(m1.bin2_ids(), m2.bin2_ids());
      compare_vectors(m1.counts(), m2.counts());
      compare_vectors(m1.chrom_offsets(), m2.chrom_offsets());
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: SparseMatrixChunked") {
  using SparseMatrixChunked = hictk::balancing::SparseMatrixChunked;
  const BinTable bins{Reference{Chromosome{0, "chr0", 50}, Chromosome{1, "chr1", 100},
                                Chromosome{2, "chr2", 50}, Chromosome{3, "chr3", 50}},
                      50};
  // clang-format off
  const std::vector<ThinPixel<std::int32_t>> pixels{
      {1, 0, 1}, {1, 1, 2}, {2, 1, 3},  // chr1
      {3, 0, 4}, {3, 1, 5}};            // chr2
  // clang-format on
  const auto tmpfile = testdir() / "sparse_matrix_chunked.tmp";

  SECTION("accessors") { CHECK(SparseMatrixChunked{bins, tmpfile, 2, 0}.empty()); }

  SECTION("push_back") {
    SparseMatrixChunked m{bins, tmpfile, 2, 0};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();

    CHECK(m.size() == pixels.size());
  }

  SECTION("subset") {
    SparseMatrixChunked m{bins, tmpfile, 2, 0};
    for (const auto& p : pixels) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
    m.finalize();

    CHECK(m.subset(0).empty());
    CHECK(m.subset(1).size() == 3);
    CHECK(m.subset(2).size() == 2);
    CHECK(m.subset(3).empty());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: ICE", "[balancing][long]") {
  const std::array<std::pair<std::string, std::filesystem::path>, 2> files{
      std::make_pair("cooler", datadir / "cooler/ENCFF993FGR.2500000.cool"),
      std::make_pair("hic", datadir / "hic/ENCFF993FGR.hic")};

  const auto tmpfile = testdir() / "balancing_ice.tmp";

  for (const auto& [label, path] : files) {
    SECTION(label) {
      const hictk::File f(path.string(), 2'500'000);

      SECTION("in-memory") {
        SECTION("INTRA") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.cis.txt";

          constexpr auto type = hictk::balancing::ICE::Type::cis;
          const auto weights = hictk::balancing::ICE(f, type).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }

        SECTION("INTER") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.trans.txt";

          constexpr auto type = hictk::balancing::ICE::Type::trans;
          const auto weights = hictk::balancing::ICE(f, type).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }

        SECTION("GW") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.gw.txt";

          constexpr auto type = hictk::balancing::ICE::Type::gw;
          const auto weights = hictk::balancing::ICE(f, type).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }
      }

      SECTION("chunked") {
        auto params = hictk::balancing::ICE::DefaultParams;
        params.tmpfile = tmpfile;
        params.chunk_size = 1000;

        SECTION("INTRA") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.cis.txt";

          constexpr auto type = hictk::balancing::ICE::Type::cis;
          const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }

        SECTION("INTER") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.trans.txt";

          constexpr auto type = hictk::balancing::ICE::Type::trans;
          const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }

        SECTION("GW") {
          const auto path_intra_weights =
              datadir / "cooler/balancing/ENCFF993FGR.2500000.ICE.gw.txt";

          constexpr auto type = hictk::balancing::ICE::Type::gw;
          const auto weights = hictk::balancing::ICE(f, type, params).get_weights();
          const auto expected_weights = read_weights(path_intra_weights);

          compare_weights(weights, expected_weights);
        }
      }
    }
  }
}

}  // namespace hictk::test::balancing
