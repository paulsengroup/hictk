// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <parallel_hashmap/btree.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/transformers/stats.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

using namespace hictk::transformers;

TEST_CASE("Transformers (cooler): stats", "[transformers][short]") {
  const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
  const cooler::File clr(path.string());
  auto sel = clr.fetch("chr1");
  auto first = sel.begin<std::int32_t>();
  auto last = sel.end<std::int32_t>();

  SECTION("range with data") {
    CHECK_THAT(avg(first, last), Catch::Matchers::WithinRel(25231.981858902574));
    CHECK(nnz(first, last) == 4'465);
    CHECK(max(first, last) == 1'357'124);
    CHECK(sum(first, last) == 112'660'799);
  }

  SECTION("empty range") {
    CHECK(avg(last, last) == 0);
    CHECK(nnz(last, last) == 0);
    CHECK(max(last, last) == 0);
    CHECK(sum(last, last) == 0);
  }
}

TEST_CASE("Transformers (hic): stats", "[transformers][short]") {
  const auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";
  const hic::File hf(path.string(), 2'500'000);
  auto sel = hf.fetch("chr2L");
  auto first = sel.begin<std::int32_t>();
  auto last = sel.end<std::int32_t>();

  SECTION("range with data") {
    CHECK_THAT(avg(first, last), Catch::Matchers::WithinRel(363057.38181818184));
    CHECK(nnz(first, last) == 55);
    CHECK(max(first, last) == 2'686'581);
    CHECK(sum(first, last) == 19'968'156);
  }

  SECTION("empty range") {
    CHECK(avg(last, last) == 0);
    CHECK(nnz(last, last) == 0);
    CHECK(max(last, last) == 0);
    CHECK(sum(last, last) == 0);
  }
}

}  // namespace hictk::test::transformers
