// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/balancing.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <vector>

#include "hictk/cooler.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::hictk {

template <typename N1, std::size_t N2>
static void balancer_test_helper(const Balancer<N1>& sel,
                                 const std::array<double, N2>& expected_counts,
                                 double abs_tol = 1.0e-6) {
  const std::vector<Pixel<double>> pixels{sel.begin(), sel.end()};

  REQUIRE(pixels.size() == expected_counts.size());
  for (std::size_t i = 0; i < expected_counts.size(); ++i) {
    CHECK_THAT(pixels[i].count, Catch::Matchers::WithinAbs(expected_counts[i], abs_tol));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: Balancer", "[cooler][short]") {
  auto path = datadir / "ENCFF993FGR.2500000.cool";
  auto clr = File::open_read_only(path.string());

  SECTION("read weights") {
    SECTION("valid") {
      CHECK(clr.read_weights("weight")->type() == Weights::Type::MULTIPLICATIVE);
      for (const auto* name : {"GW_SCALE", "INTER_SCALE", "SCALE", "VC", "VC_SQRT"}) {
        CHECK(clr.read_weights(name)->type() == Weights::Type::DIVISIVE);
      }
    }

    SECTION("invalid") {
      CHECK_THROWS(clr.read_weights(""));
      CHECK_THROWS(clr.read_weights("AAA"));
    }

    SECTION("purging") {
      CHECK(clr.purge_weights() == false);
      CHECK(clr.purge_weights("weight") == false);

      const auto w = clr.read_weights("weight");
      CHECK(w.use_count() == 2);
      CHECK(clr.purge_weights("weight") == true);
      CHECK(w.use_count() == 1);

      clr.read_weights("weight");
      CHECK(clr.purge_weights() == true);
    }
  }

  SECTION("Balancer") {
    SECTION("operators") {
      const auto sel = Balancer(clr.fetch<std::int32_t>("chr1", 5'000'000, 10'000'000),
                                clr.read_weights("weight"));

      CHECK(sel.begin() < ++sel.begin());
      CHECK(sel.begin() <= ++sel.begin());
      CHECK(sel.begin() <= sel.begin()++);
      CHECK(++sel.begin() > sel.begin());
      CHECK(++sel.begin() >= sel.begin());
      CHECK(sel.begin() >= sel.begin()++);

      CHECK_THAT(sel.begin()->count, Catch::Matchers::WithinAbs(3.345797, 1.0e-6));
    }

    SECTION("cis") {
      SECTION("ICE") {
        const auto sel = Balancer(clr.fetch<std::int32_t>("chr1", 5'000'000, 10'000'000),
                                  clr.read_weights("weight"));
        constexpr std::array<double, 3> expected{3.345797, 0.328794, 4.456354};
        balancer_test_helper(sel, expected);
      }

      SECTION("GW_SCALE") {
        const auto sel = Balancer(clr.fetch<std::int32_t>("chr1", 5'000'000, 10'000'000),
                                  clr.read_weights("GW_SCALE"));
        constexpr std::array<double, 3> expected{927703.336647, 77376.912375, 890112.397104};
        balancer_test_helper(sel, expected);
      }
    }

    SECTION("trans") {
      SECTION("ICE") {
        const auto sel = Balancer(
            clr.fetch<std::int32_t>("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000),
            clr.read_weights("weight"));
        constexpr std::array<double, 4> expected{0.001782, 0.002756, 0.002047, 0.004749};
        balancer_test_helper(sel, expected);
      }

      SECTION("GW_SCALE") {
        const auto sel = Balancer(
            clr.fetch<std::int32_t>("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000),
            clr.read_weights("GW_SCALE"));
        constexpr std::array<double, 4> expected{600.616151, 761.596365, 585.635384, 1113.900564};
        balancer_test_helper(sel, expected);
      }
    }
  }
}
}  // namespace hictk::test::hictk
