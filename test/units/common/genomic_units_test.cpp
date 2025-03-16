// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/genomic_units.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

namespace hictk::test::common {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

template <typename N>
[[nodiscard]] static std::string format_distance(N n, std::string_view multiplier = "",
                                                 std::string_view suffix = "",
                                                 bool with_space = false) {
  return fmt::format(FMT_STRING("{}{}{}{}"), n, with_space ? " " : "", multiplier, suffix);
}

TEST_CASE("Common: parse_genomic_unit", "[common][short]") {
  SECTION("valid") {
    CHECK(parse_genomic_unit("bp") == 1);
    CHECK(parse_genomic_unit("kbp") == 1'000);
    CHECK(parse_genomic_unit("mbp") == 1'000'000);
    CHECK(parse_genomic_unit("gbp") == 1'000'000'000);
  }

  SECTION("invalid") {
    CHECK_THROWS_WITH(parse_genomic_unit(""), Catch::Matchers::ContainsSubstring("is empty"));
    CHECK_THROWS_WITH(parse_genomic_unit("abcd"),
                      Catch::Matchers::ContainsSubstring("Unrecognized unit"));
    CHECK_THROWS_WITH(parse_genomic_unit("kx"),
                      Catch::Matchers::ContainsSubstring("Unrecognized unit"));
    CHECK_THROWS_WITH(parse_genomic_unit("kxp"),
                      Catch::Matchers::ContainsSubstring("Unrecognized unit"));
    CHECK_THROWS_WITH(parse_genomic_unit("kbx"),
                      Catch::Matchers::ContainsSubstring("Unrecognized unit"));
    CHECK_THROWS_WITH(parse_genomic_unit("tbp"),
                      Catch::Matchers::ContainsSubstring("Unrecognized unit"));
  }
}

TEST_CASE("Common: parse_genomic_distance", "[common][short]") {
  SECTION("valid") {
    for (const std::string_view suffix : {"", "b", "bp", "B", "BP"}) {
      if (suffix != "b" && suffix != "B") {
        CHECK(parse_genomic_distance(format_distance(2)) == 2);
        CHECK(parse_genomic_distance(format_distance(2, "", suffix)) == 2);
      }

      CHECK(parse_genomic_distance(format_distance(2, "k", suffix)) == 2'000);
      CHECK(parse_genomic_distance(format_distance(2, "K", suffix)) == 2'000);

      CHECK(parse_genomic_distance(format_distance(2, "m", suffix)) == 2'000'000);
      CHECK(parse_genomic_distance(format_distance(2, "M", suffix)) == 2'000'000);

      CHECK(parse_genomic_distance(format_distance(2, "g", suffix)) == 2'000'000'000);
      CHECK(parse_genomic_distance(format_distance(2, "G", suffix)) == 2'000'000'000);

      CHECK(parse_genomic_distance(format_distance(2.1, "k", suffix)) == 2'100);
      CHECK(parse_genomic_distance(format_distance(2.1, "K", suffix)) == 2'100);

      CHECK(parse_genomic_distance(format_distance(2.1, "m", suffix)) == 2'100'000);
      CHECK(parse_genomic_distance(format_distance(2.1, "M", suffix)) == 2'100'000);

      CHECK(parse_genomic_distance(format_distance(2.1, "g", suffix)) == 2'100'000'000);
      CHECK(parse_genomic_distance(format_distance(2.1, "G", suffix)) == 2'100'000'000);
    }
  }

  SECTION("invalid") {
    CHECK_THROWS_WITH(parse_genomic_distance(""), Catch::Matchers::ContainsSubstring("is empty"));
    CHECK_THROWS_WITH(parse_genomic_distance(".123"),
                      Catch::Matchers::ContainsSubstring("does not start with a digit"));
    CHECK_THROWS_WITH(parse_genomic_distance("a123"),
                      Catch::Matchers::ContainsSubstring("does not start with a digit"));
    CHECK_THROWS(parse_genomic_distance("123.123.123"));
    CHECK_THROWS_WITH(parse_genomic_distance("123 "),
                      Catch::Matchers::ContainsSubstring("has trailing whitespaces"));
    CHECK_THROWS_WITH(parse_genomic_distance("1.2345 kbp"),
                      Catch::Matchers::ContainsSubstring("Cannot convert") &&
                          Catch::Matchers::ContainsSubstring("to an integer"));
    CHECK_THROWS(parse_genomic_distance<std::int8_t>("200"));
    CHECK_THROWS_WITH(parse_genomic_distance<std::int8_t>("1 kbp"),
                      Catch::Matchers::ContainsSubstring("Cannot fit"));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::common
