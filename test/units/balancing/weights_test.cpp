// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/weights.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <numeric>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

namespace hictk::test::balancing {

using Weights = hictk::balancing::Weights;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weights (vector)", "[balancing][weights][short]") {
  const std::vector<double> weights{1.0, 2.0, 3.0, 4.0, 5.0};

  SECTION("Ctors") {
    auto w = Weights(weights, Weights::Type::MULTIPLICATIVE);
    CHECK(w.type() == Weights::Type::MULTIPLICATIVE);
    CHECK(w.size() == weights.size());

    CHECK_THROWS(Weights(weights, Weights::Type::UNKNOWN));

    for (const auto& name : {"VC", "INTER_VC", "GW_VC", "VC_SQRT", "KR", "INTER_KR", "GW_KR",
                             "SCALE", "INTER_SCALE", "GW_SCALE"}) {
      w = Weights(weights, name);
      CHECK(w.type() == Weights::Type::DIVISIVE);
      CHECK(w.size() == weights.size());
    }

    for (const auto& name : {"ICE", "INTER_ICE", "GW_ICE", "weight"}) {
      w = Weights(weights, name);
      CHECK(w.type() == Weights::Type::MULTIPLICATIVE);
      CHECK(w.size() == weights.size());
    }

    CHECK_THROWS(Weights(weights, "foobar"));
  }

  const Weights wm(weights, Weights::Type::MULTIPLICATIVE);
  const Weights wd(weights, Weights::Type::DIVISIVE);

  SECTION("accessors") {
    SECTION("operator bool()") {
      CHECK(!Weights{});
      CHECK(!!wm);
      CHECK(!!wd);
    }

    SECTION("operator[]") {
      REQUIRE(wm.size() == weights.size());
      REQUIRE(wd.size() == weights.size());
      for (std::size_t i = 0; i < weights.size(); ++i) {
        CHECK(wm[i] == weights[i]);
        CHECK(wd[i] == weights[i]);
      }
    }

    SECTION("at()") {
      CHECK_THROWS(Weights{}.at(0));
      CHECK_THROWS(Weights{}.at(5));

      for (std::size_t i = 0; i < weights.size(); ++i) {
        CHECK(wm.at(i) == weights.at(i));
        CHECK(wd.at(i) == weights.at(i));

        CHECK(wm.at(i, Weights::Type::MULTIPLICATIVE) == weights.at(i));
        CHECK(wd.at(i, Weights::Type::DIVISIVE) == weights.at(i));

        CHECK(wm.at(i, Weights::Type::DIVISIVE) == 1.0 / weights.at(i));
        CHECK(wd.at(i, Weights::Type::MULTIPLICATIVE) == 1.0 / weights.at(i));
      }

      CHECK_THROWS(wm.at(weights.size() + 1));
      CHECK_THROWS(wm.at(0, Weights::Type::UNKNOWN));
    }

    SECTION("to_vector()") {
      auto v = wm.to_vector();
      CHECK(std::accumulate(v.begin(), v.end(), 0.0) == std::accumulate(wm.begin(), wm.end(), 0.0));
      v = wm.to_vector(Weights::Type::DIVISIVE);
      CHECK(
          std::accumulate(v.begin(), v.end(), 0.0) ==
          std::accumulate(wm.begin(Weights::Type::DIVISIVE), wm.end(Weights::Type::DIVISIVE), 0.0));
    }
  }

  SECTION("balance()") {
    CHECK(wm.balance(ThinPixel<double>{0, 0, 1}).count == 1.0);
    CHECK(wm.balance(ThinPixel<double>{3, 4, 1}).count == 20.0);

    CHECK(wd.balance(ThinPixel<double>{0, 0, 1}).count == 1.0);
    CHECK(wd.balance(ThinPixel<double>{3, 4, 1}).count == 1.0 / 20.0);
  }

  SECTION("operator()") {
    auto w = wm(Weights::Type::MULTIPLICATIVE);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == wm.at(i));
    }

    w = wm(Weights::Type::DIVISIVE);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == 1.0 / wm.at(i));
    }

    w = wd(Weights::Type::DIVISIVE);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == wd.at(i));
    }

    w = wd(Weights::Type::MULTIPLICATIVE);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == 1.0 / wd.at(i));
    }

    CHECK_THROWS(wm(Weights::Type::INFER));
  }

  SECTION("rescale()") {
    auto w = wm;
    w.rescale(2.0);
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == std::sqrt(2.0) * wm.at(i));
    }

    w = wm;
    w.rescale({2.0}, {0, weights.size()});
    for (std::size_t i = 0; i < weights.size(); ++i) {
      CHECK(w.at(i) == std::sqrt(2.0) * wm.at(i));
    }

    const std::vector<double> sclaing_factors{1.0, 10.0};
    const std::vector<std::uint64_t> offsets{0, 3, weights.size()};

    w = wm;
    w.rescale(sclaing_factors, offsets);

    CHECK(w.at(0) == wm.at(0));
    CHECK(w.at(1) == wm.at(1));
    CHECK(w.at(2) == wm.at(2));
    CHECK(w.at(3) == std::sqrt(10.0) * wm.at(3));
    CHECK(w.at(4) == std::sqrt(10.0) * wm.at(4));

    CHECK_THROWS(w.rescale({}, {}));
    CHECK_THROWS(w.rescale({1}, {0}));                              // invalid offsets size
    CHECK_THROWS(w.rescale({1}, {1, weights.size()}));              // invalid offsets.front()
    CHECK_THROWS(w.rescale({1}, {0, 1}));                           // invalid offsets.back()
    CHECK_THROWS(w.rescale({1, 1, 1}, {0, 2, 1, weights.size()}));  // offsets not sorted
  }

  SECTION("iteration") {
    std::for_each(wm.begin(), wm.end(),
                  [&, i = std::size_t{}](const auto& w) mutable { CHECK(w == wm.at(i++)); });
    std::for_each(wm.begin(Weights::Type::DIVISIVE), wm.end(Weights::Type::DIVISIVE),
                  [&, i = std::size_t{}](const auto& w) mutable { CHECK(1.0 / w == wm.at(i++)); });

    const auto it1 = wm.begin();
    const auto it2 = wm.begin() + 1;
    const auto it3 = wm.begin() + 2;

    CHECK(it1 == it1);
    CHECK(it1 != it2);
    CHECK(it1 < it2);
    CHECK(it1 <= it1);
    CHECK(it1 <= it2);
    CHECK_FALSE(it1 > it2);
    CHECK_FALSE(it1 >= it2);
    CHECK(it2 >= it2);

    CHECK(it1[2] == it2[1]);
    CHECK(it1[0] == it3[-2]);

    CHECK(++wm.begin() == it2);
    CHECK(wm.begin()++ == it1);

    auto it = it1;
    it += 2;
    CHECK(it == it3);

    CHECK(it1 + 1 == it2);

    it = it2;
    CHECK(it1 == --it);

    it = it2;
    it -= 1;
    CHECK(it1 == it);

    CHECK(it1 == it2 - 1);

    CHECK(it3 - it1 == 2);

    CHECK_THROWS(wm.begin(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.end(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.cbegin(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.cend(Weights::Type::UNKNOWN));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Weights (constant)", "[balancing][weights][short]") {
  SECTION("Ctors") {
    auto w = Weights(1.0, 5, Weights::Type::MULTIPLICATIVE);
    CHECK(w.type() == Weights::Type::MULTIPLICATIVE);
    CHECK(w.size() == 5);

    CHECK_THROWS(Weights(1.0, 5.0, Weights::Type::UNKNOWN));

    for (const auto& name : {"VC", "INTER_VC", "GW_VC", "VC_SQRT", "KR", "INTER_KR", "GW_KR",
                             "SCALE", "INTER_SCALE", "GW_SCALE"}) {
      w = Weights(1.0, 5, name);
      CHECK(w.type() == Weights::Type::DIVISIVE);
      CHECK(w.size() == 5);
    }

    for (const auto& name : {"ICE", "INTER_ICE", "GW_ICE", "weight"}) {
      w = Weights(1.0, 5, name);
      CHECK(w.type() == Weights::Type::MULTIPLICATIVE);
      CHECK(w.size() == 5);
    }

    CHECK_THROWS(Weights(1.0, 5, "foobar"));
  }

  const Weights wm(5.0, 10, Weights::Type::MULTIPLICATIVE);
  const Weights wd(5.0, 10, Weights::Type::DIVISIVE);

  SECTION("accessors") {
    SECTION("operator bool()") {
      CHECK(!Weights{});
      CHECK(!!wm);
      CHECK(!!wd);
    }

    SECTION("operator[]") {
      REQUIRE(wm.size() == 10);
      REQUIRE(wd.size() == 10);
      for (std::size_t i = 0; i < 10; ++i) {
        CHECK(wm[i] == 5.0);
        CHECK(wd[i] == 5.0);
      }
    }

    SECTION("at()") {
      CHECK_THROWS(Weights{}.at(0));
      CHECK_THROWS(Weights{}.at(5));

      for (std::size_t i = 0; i < 10; ++i) {
        CHECK(wm.at(i) == 5.0);
        CHECK(wd.at(i) == 5.0);

        CHECK(wm.at(i, Weights::Type::MULTIPLICATIVE) == 5.0);
        CHECK(wd.at(i, Weights::Type::DIVISIVE) == 5.0);

        CHECK(wm.at(i, Weights::Type::DIVISIVE) == 1.0 / 5.0);
        CHECK(wd.at(i, Weights::Type::MULTIPLICATIVE) == 1.0 / 5.0);
      }

      CHECK_THROWS(wm.at(11));
      CHECK_THROWS(wm.at(0, Weights::Type::UNKNOWN));
    }

    SECTION("to_vector()") {
      auto v = wm.to_vector();
      CHECK(std::accumulate(v.begin(), v.end(), 0.0) == std::accumulate(wm.begin(), wm.end(), 0.0));
      v = wm.to_vector(Weights::Type::DIVISIVE);
      CHECK(
          std::accumulate(v.begin(), v.end(), 0.0) ==
          std::accumulate(wm.begin(Weights::Type::DIVISIVE), wm.end(Weights::Type::DIVISIVE), 0.0));
    }
  }

  SECTION("balance()") {
    CHECK(wm.balance(ThinPixel<double>{0, 0, 1}).count == 25.0);
    CHECK(wm.balance(ThinPixel<double>{3, 4, 1}).count == 25.0);

    CHECK(wd.balance(ThinPixel<double>{0, 0, 1}).count == 1.0 / 25.0);
    CHECK(wd.balance(ThinPixel<double>{3, 4, 1}).count == 1.0 / 25.0);
  }

  SECTION("operator()") {
    auto w = wm(Weights::Type::MULTIPLICATIVE);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == wm.at(i));
    }

    w = wm(Weights::Type::DIVISIVE);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == 1.0 / wm.at(i));
    }

    w = wd(Weights::Type::DIVISIVE);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == wd.at(i));
    }

    w = wd(Weights::Type::MULTIPLICATIVE);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == 1.0 / wd.at(i));
    }

    CHECK_THROWS(wm(Weights::Type::INFER));
  }

  SECTION("rescale()") {
    auto w = wm;
    w.rescale(2.0);
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == std::sqrt(2.0) * wm.at(i));
    }

    w = wm;
    w.rescale({2.0}, {0, 10});
    for (std::size_t i = 0; i < 10; ++i) {
      CHECK(w.at(i) == std::sqrt(2.0) * wm.at(i));
    }

    CHECK_THROWS(w.rescale({}, {}));
    CHECK_THROWS(w.rescale({1}, {0}));            // invalid offsets size
    CHECK_THROWS(w.rescale({1}, {1, 10}));        // invalid offsets.front()
    CHECK_THROWS(w.rescale({1}, {0, 1}));         // invalid offsets.back()
    CHECK_THROWS(w.rescale({1, 2}, {0, 5, 10}));  // invalid shape
  }

  SECTION("iteration") {
    std::for_each(wm.begin(), wm.end(),
                  [&, i = std::size_t{}](const auto& w) mutable { CHECK(w == wm.at(i++)); });
    std::for_each(wm.begin(Weights::Type::DIVISIVE), wm.end(Weights::Type::DIVISIVE),
                  [&, i = std::size_t{}](const auto& w) mutable { CHECK(1.0 / w == wm.at(i++)); });

    const auto it1 = wm.begin();
    const auto it2 = wm.begin() + 1;
    const auto it3 = wm.begin() + 2;

    CHECK(it1 == it1);
    CHECK(it1 != it2);
    CHECK(it1 < it2);
    CHECK(it1 <= it1);
    CHECK(it1 <= it2);
    CHECK_FALSE(it1 > it2);
    CHECK_FALSE(it1 >= it2);
    CHECK(it2 >= it2);

    CHECK(it1[5] == it2[4]);
    CHECK(it1[0] == it2[-1]);

    CHECK(++wm.begin() == it2);
    CHECK(wm.begin()++ == it1);

    auto it = it1;
    it += 2;
    CHECK(it == it3);

    CHECK(it1 + 1 == it2);

    it = it2;
    CHECK(it1 == --it);

    it = it2;
    it -= 1;
    CHECK(it1 == it);

    CHECK(it1 == it2 - 1);

    CHECK(it3 - it1 == 2);

    CHECK_THROWS(wm.begin(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.end(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.cbegin(Weights::Type::UNKNOWN));
    CHECK_THROWS(wm.cend(Weights::Type::UNKNOWN));
  }
}

}  // namespace hictk::test::balancing
