// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <random>
#include <variant>
#include <vector>

#include "hictk/variant_buff.hpp"

namespace hictk::test::variantbuff {
using namespace hictk::internal;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("VariantBuffer: ctor", "[numeric_buff][short]") {
  const std::vector<double> buff0(10, 1.0);

  SECTION("default ctor") {
    const VariantBuffer buff1{};

    CHECK(buff1.empty());
  }

  SECTION("ctor 1") {
    VariantBuffer buff1(10, 1.0);

    CHECK(buff0.size() == buff1.size<double>());
    CHECK(buff0.size() == buff1.size());

    CHECK(buff0.capacity() == buff1.capacity<double>());
    CHECK(buff0.capacity() == buff1.capacity());

    CHECK(std::equal(buff1.begin<double>(), buff1.end<double>(), buff0.begin()));

    buff1.clear<double>();
    CHECK(buff1.empty<double>());
    CHECK(buff1.empty());
  }

  SECTION("ctor 2") {
    VariantBuffer buff1(buff0);

    CHECK(buff0.size() == buff1.size<double>());
    CHECK(buff0.size() == buff1.size());

    CHECK(buff0.capacity() == buff1.capacity<double>());
    CHECK(buff0.capacity() == buff1.capacity());

    CHECK(std::equal(buff1.begin<double>(), buff1.end<double>(), buff0.begin()));

    buff1.clear<double>();
    CHECK(buff1.empty<double>());
    CHECK(buff1.empty());
  }

  SECTION("ctor 3") {
    VariantBuffer buff1(buff0.begin(), buff0.end());

    CHECK(buff0.size() == buff1.size<double>());
    CHECK(buff0.size() == buff1.size());

    CHECK(buff0.capacity() == buff1.capacity<double>());
    CHECK(buff0.capacity() == buff1.capacity());

    CHECK(std::equal(buff1.begin<double>(), buff1.end<double>(), buff0.begin()));

    buff1.clear<double>();
    CHECK(buff1.empty<double>());
    CHECK(buff1.empty());
  }
}

TEST_CASE("VariantBuffer: accessors", "[numeric_buff][short]") {
  using T = std::uint64_t;

  std::vector<T> buff0(10);
  std::iota(buff0.begin(), buff0.end(), 0);

  const VariantBuffer buff1(buff0);

  SECTION("bad variant access") { CHECK_THROWS_AS(buff1.get<int>(), std::bad_variant_access); }

  SECTION("front(), back() and data()") {
    CHECK(buff0.front() == buff1.front<T>());
    CHECK(buff0.back() == buff1.back<T>());
    CHECK(*buff0.data() == *buff1.data<T>());
  }

  SECTION("at() and operator[]") {
    for (std::size_t i = 0; i < buff0.size(); ++i) {
      CHECK(buff0[i] == buff1.at<T>(i));
      CHECK(buff0[i] == std::get<T>(buff1.at(i)));
      CHECK(buff0[i] == std::get<T>(buff1[i]));
    }
  }

  SECTION("begin()/end()") {
    std::mt19937_64 prng{};  // NOLINT
    auto buff2 = buff1;

    std::shuffle(buff2.begin<std::uint64_t>(), buff2.end<std::uint64_t>(), prng);

    CHECK(
        std::accumulate(buff0.begin(), buff0.end(), std::size_t(0)) ==
        std::accumulate(buff2.begin<std::uint64_t>(), buff2.end<std::uint64_t>(), std::size_t(0)));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::variantbuff
