// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/pixel.hpp"

#include <catch2/catch_test_macros.hpp>

#include "coolerpp/chromosome.hpp"

namespace coolerpp {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Pixel", "[pixel][short]") {
  const ChromosomeSet chroms{
      Chromosome{0, "chr1", 248'956'422}, Chromosome{1, "chr2", 242'193'529},
      Chromosome{2, "chr3", 198'295'559}, Chromosome{3, "chr4", 190'214'555},
      Chromosome{4, "chr5", 181'538'259}, Chromosome{5, "chr6", 170'805'979},
      Chromosome{6, "chr9", 138'394'717}, Chromosome{7, "chr11", 135'086'622},
      Chromosome{8, "chr12", 133'275'309}};
  constexpr std::uint32_t bin_size = 1;
  const auto bins = std::make_shared<const BinTable>(chroms, bin_size);

  auto PI = [&](std::string_view chrom1, std::string_view chrom2, std::uint32_t pos1,
                std::uint32_t pos2, std::int32_t count = 0) {
    return Pixel<std::int32_t>{Bin{chroms.at(chrom1), pos1, pos1 + bin_size},
                               Bin{chroms.at(chrom2), pos2, pos2 + bin_size}, count};
  };

  auto PFP = [&](std::string_view chrom1, std::string_view chrom2, std::uint32_t pos1,
                 std::uint32_t pos2, double count = 0) {
    return Pixel<double>{Bin{chroms.at(chrom1), pos1, pos1 + bin_size},
                         Bin{chroms.at(chrom2), pos2, pos2 + bin_size}, count};
  };

  SECTION("operator bool") {
    CHECK(!Pixel<std::int32_t>{});
    CHECK(!!PI("chr1", "chr1", 0, 10));
  }

  SECTION("(dis)equality") {
    CHECK(PI("chr1", "chr1", 0, 10) == PI("chr1", "chr1", 0, 10));

    CHECK(PI("chr1", "chr1", 0, 10) != PI("chr1", "chr2", 0, 10));
    CHECK(PI("chr1", "chr1", 0, 10) != PI("chr2", "chr1", 0, 10));

    CHECK(PI("chr1", "chr1", 0, 10) != PI("chr1", "chr1", 1, 10));
    CHECK(PI("chr1", "chr1", 0, 10) != PI("chr1", "chr1", 0, 11));
  }

  SECTION("ordering") {
    CHECK(PI("chr1", "chr1", 0, 0) < PI("chr2", "chr2", 0, 0));
    CHECK(PI("chr1", "chr1", 0, 0) <= PI("chr2", "chr2", 0, 0));

    CHECK(PI("chr1", "chr1", 0, 0) < PI("chr1", "chr2", 0, 0));
    CHECK(PI("chr1", "chr1", 0, 0) <= PI("chr1", "chr2", 0, 0));

    CHECK(PI("chr2", "chr2", 0, 0) > PI("chr1", "chr1", 0, 0));
    CHECK(PI("chr2", "chr2", 0, 0) >= PI("chr1", "chr1", 0, 0));

    CHECK(PI("chr1", "chr2", 0, 0) > PI("chr1", "chr1", 0, 0));
    CHECK(PI("chr1", "chr2", 0, 0) >= PI("chr1", "chr1", 0, 0));

    CHECK(PI("chr1", "chr1", 0, 0) < PI("chr1", "chr1", 0, 1));
    CHECK(PI("chr1", "chr1", 0, 0) < PI("chr1", "chr1", 1, 0));
    CHECK(PI("chr1", "chr1", 0, 0) <= PI("chr1", "chr1", 0, 1));
    CHECK(PI("chr1", "chr1", 0, 0) <= PI("chr1", "chr1", 1, 0));

    CHECK(PI("chr1", "chr1", 0, 1) > PI("chr1", "chr1", 0, 0));
    CHECK(PI("chr1", "chr1", 1, 0) > PI("chr1", "chr1", 0, 0));
    CHECK(PI("chr1", "chr1", 0, 1) >= PI("chr1", "chr1", 0, 0));
    CHECK(PI("chr1", "chr1", 1, 0) >= PI("chr1", "chr1", 0, 0));
  }

  SECTION("sorting") {
    const std::vector<Pixel<std::int32_t>> pixels{{
        // clang-format off
        PI("chr1", "chr1",  10'000, 180'000),
        PI("chr1", "chr1",  10'000, 202'890'000),
        PI("chr1", "chr2",  10'000, 113'590'000),
        PI("chr1", "chr4",  10'000, 52'880'000),
        PI("chr1", "chr5",  10'000, 230'000),
        PI("chr1", "chr6",  10'000, 33'820'000),
        PI("chr1", "chr6",  10'000, 149'280'000),
        PI("chr1", "chr9",  10'000, 10'000),
        PI("chr1", "chr9",  10'000, 122'380'000),
        PI("chr1", "chr11", 40'000, 11'630'000),
        PI("chr1", "chr11", 40'000, 120'770'000),
        PI("chr1", "chr12", 40'000, 7'060'000),
        PI("chr1", "chr12", 40'000, 119'750'000),
        PI("chr2", "chr2",  10'000, 10'000),
        PI("chr2", "chr2",  10'000, 20'000),
        PI("chr2", "chr3",  10'000, 99'320'000),
        PI("chr2", "chr3",  10'000, 101'660'000),
        // clang-format on
    }};

    CHECK(std::is_sorted(pixels.begin(), pixels.end()));
  }

  SECTION("fmt") {
    auto p1 = PI("chr1", "chr1", 0, 10);
    REQUIRE(p1.coords.bin1.has_null_id());
    REQUIRE(p1.coords.bin2.has_null_id());

    CHECK(fmt::format(FMT_STRING("{}"), p1) == "chr1\t0\t1\tchr1\t10\t11\t0");
    CHECK(fmt::format(FMT_STRING("{:bg2}"), p1) == "chr1\t0\t1\tchr1\t10\t11\t0");
    CHECK(fmt::format(FMT_STRING("{:raw}"), p1) == "18446744073709551615\t18446744073709551615\t0");

    auto p2 = PFP("chr1", "chr1", 0, 10, 1.2);
    REQUIRE(p2.coords.bin1.has_null_id());
    REQUIRE(p2.coords.bin2.has_null_id());
    CHECK(fmt::format(FMT_STRING("{}"), p2) == "chr1\t0\t1\tchr1\t10\t11\t1.2");
    CHECK(fmt::format(FMT_STRING("{:bg2}"), p2) == "chr1\t0\t1\tchr1\t10\t11\t1.2");
    CHECK(fmt::format(FMT_STRING("{:raw}"), p2) ==
          "18446744073709551615\t18446744073709551615\t1.2");
  }
}

}  // namespace coolerpp
