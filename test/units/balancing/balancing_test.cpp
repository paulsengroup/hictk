// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <zstd.h>

#include <array>
#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "./common.hpp"
#include "hictk/balancing/ice.hpp"
#include "hictk/balancing/scale.hpp"
#include "hictk/balancing/vc.hpp"
#include "hictk/file.hpp"
#include "tmpdir.hpp"

namespace hictk::test {
// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

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
        std::filesystem::remove(tmpfile);  // NOLINT
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
        std::filesystem::remove(tmpfile);  // NOLINT
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
        std::filesystem::remove(tmpfile);  // NOLINT
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
        std::filesystem::remove(tmpfile);  // NOLINT
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
        std::filesystem::remove(tmpfile);  // NOLINT
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
        std::filesystem::remove(tmpfile);  // NOLINT
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
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::balancing
