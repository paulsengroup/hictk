// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>
#include <set>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset read", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      constexpr std::array<std::string_view, 3> expected{"1", "2", "3"};
      std::vector<std::string> buff{};

      Dataset{grp, "chroms/name"}.read(buff, expected.size());
      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/name"};
      std::string buff;
      dset.read(buff, 9);
      CHECK(buff == "10");
      CHECK(dset.read_last<std::string>() == "X");
      CHECK(std::get<std::string>(dset.read_last()) == "X");
    }
  }

  SECTION("numeric") {
    using T = std::int32_t;
    constexpr std::array<T, 10> expected{0,       100'000, 200'000, 300'000, 400'000,
                                         500'000, 600'000, 700'000, 800'000, 900'000};

    constexpr std::size_t nnz_expected = 107'041;
    constexpr std::int32_t sum_expected = 395'465;

    SECTION("vector<T>") {
      std::vector<T> buff{};
      std::ignore = Dataset{grp, "bins/start"}.read(buff, expected.size());

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(buff);
      CHECK(buff.size() == nnz_expected);
      CHECK(std::accumulate(buff.begin(), buff.end(), 0) == sum_expected);
    }

    SECTION("variant buff") {
      hictk::internal::VariantBuffer vbuff{std::size_t(0), 0.0};
      std::ignore = Dataset{grp, "bins/start"}.read(vbuff, expected.size());
      const auto& buff = vbuff.get<T>();
      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(vbuff);
      CHECK(vbuff.size<T>() == nnz_expected);
      CHECK(std::accumulate(vbuff.begin<T>(), vbuff.end<T>(), 0) == sum_expected);
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/length"};
      std::uint64_t buff{};
      dset.read(buff, 2);
      CHECK(buff == 159'599'783);

      CHECK(dset.read_last<std::int32_t>() == 166'650'296);
      CHECK(std::get<std::int32_t>(dset.read_last()) == 166'650'296);
    }

    SECTION("enum") { CHECK(Dataset{grp, "bins/chrom"}.read<std::uint32_t>(0) == 0); }
  }
}

}  // namespace hictk::cooler::test::dataset
