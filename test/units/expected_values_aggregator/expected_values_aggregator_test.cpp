// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/expected_values_aggregator.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/hic.hpp"

using namespace hictk;

namespace hictk::test::expected_values_aggregator {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("ExpectedValuesAggregator", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();

  const hictk::hic::File f(path_hic, resolution);

  ExpectedValuesAggregator aggr(f.bins_ptr());

  const auto sel = f.fetch();
  std::for_each(sel.template begin<std::uint32_t>(), sel.template end<std::uint32_t>(),
                [&](const auto& tp) { aggr.add(tp); });
  aggr.compute_density();

  SECTION("valid chromosome") {
    const auto chrom = f.chromosomes().longest_chromosome();

    const auto expected = f.expected_values(chrom);
    const auto expected_computed = aggr.weights(chrom);

    REQUIRE(expected.size() == expected_computed.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK_THAT(expected_computed[i], Catch::Matchers::WithinRel(expected[i]));
    }
  }

  SECTION("invalid chromosome") { CHECK_THROWS(aggr.weights(Chromosome{99, "A", 10})); }

  SECTION("small chromosome") {
    const BinTable bins{{Chromosome{0, "chr1", 5}}, 10};

    ExpectedValuesAggregator aggr1(std::make_shared<const BinTable>(bins));
    CHECK_NOTHROW(aggr1.compute_density());
    CHECK(aggr1.weights().empty());
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::expected_values_aggregator
