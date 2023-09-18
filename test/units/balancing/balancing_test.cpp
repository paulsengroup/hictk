// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/vc.hpp"
#include "hictk/hic.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::balancing {

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Balancing: VC", "[balancing][short]") {
  const auto path = datadir / "ENCFF993FGR.hic";

  auto hf = hictk::hic::File(path.string(), 2500000);

  SECTION("INTRA") {
    for (const auto& chrom : hf.chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      auto sel1 = hf.fetch(chrom.name());

      const auto num_bins = hf.bins().subset(chrom).size();
      const auto bin_id_offset = hf.bins().at(chrom.name(), 0).id();
      const auto weights =
          hictk::balancing::VC<std::int32_t>(sel1.begin<std::int32_t>(), sel1.end<std::int32_t>(),
                                             num_bins, bin_id_offset)
              .get_weights();

      auto sel2 = hf.fetch(chrom.name(), hictk::balancing::Method::VC());
      compare_weights(weights, sel2.weights1()());
    }
  }

  SECTION("GW") {
    const auto num_bins = hf.bins().size();
    auto sel = hf.fetch();
    const auto weights = hictk::balancing::VC<std::int32_t>(sel.begin<std::int32_t>(),
                                                            sel.end<std::int32_t>(), num_bins)
                             .get_weights();

    const auto expected = hf.fetch(hictk::balancing::Method::GW_VC()).weights();
    compare_weights(weights, expected);
  }
}

}  // namespace hictk::test::balancing
