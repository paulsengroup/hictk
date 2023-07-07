// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/attribute.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <highfive/H5File.hpp>
#include <string>

#include "hictk/suppress_warnings.hpp"
#include "tmpdir.hpp"


namespace hictk::cooler::test::attribute {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: attribute read - test numeric conversions", "[cooler][short]") {
  const auto path = testdir() / "test_read_attrs_numeric_conversion.cool";

  auto f = HighFive::File(path.string(), HighFive::File::Truncate);

  const double dbl = 10.0;
  const float flt = 10.0;
  const std::int32_t i32 = 12345;
  const std::string dbl_str = std::to_string(dbl);
  const std::string int_str = std::to_string(i32);

  Attribute::write(f, "double", dbl);
  Attribute::write(f, "float", flt);
  Attribute::write(f, "std::int32_t", i32);
  Attribute::write(f, "double_s", dbl_str);
  Attribute::write(f, "int_s", int_str);

  SECTION("no conversion") { CHECK(Attribute::read<double>(f, "double") == dbl); }

  SECTION("double to float") { CHECK(Attribute::read<float>(f, "double") == Catch::Approx(flt)); }
  SECTION("float to double") { CHECK(Attribute::read<double>(f, "float") == Catch::Approx(dbl)); }

  SECTION("int lossless") {
    CHECK(Attribute::read<std::uint32_t>(f, "std::int32_t") == i32);

    Attribute::write(f, "std::int64_t", std::int64_t(-1));
    CHECK(Attribute::read<std::int8_t>(f, "std::int64_t") == -1);
  }

  SECTION("int lossy") {
    CHECK_THROWS_WITH(Attribute::read<std::int8_t>(f, "std::int32_t"),
                      Catch::Matchers::ContainsSubstring("Unable to represent value 12345") &&
                          Catch::Matchers::ContainsSubstring("without overflowing"));

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_USELESS_CAST
    Attribute::write(f, "std::int32_t", std::int32_t(-1), true);
    CHECK_THROWS_WITH(Attribute::read<std::uint64_t>(f, "std::int32_t"),
                      Catch::Matchers::ContainsSubstring("Unable to represent value -1") &&
                          Catch::Matchers::ContainsSubstring("without overflowing"));
    DISABLE_WARNING_POP

    Attribute::write(f, "std::int64_t", std::numeric_limits<std::int64_t>::lowest());
    CHECK_THROWS_WITH(Attribute::read<std::int32_t>(f, "std::int64_t"),
                      Catch::Matchers::ContainsSubstring("Unable to represent value") &&
                          Catch::Matchers::ContainsSubstring("without overflowing"));
  }

  SECTION("str to double") { CHECK(Attribute::read<double>(f, "double_s") == dbl); }

  SECTION("str to int lossless") { CHECK(Attribute::read<std::int32_t>(f, "int_s") == i32); }

  SECTION("str to int lossy") {
    CHECK_THROWS_WITH(
        Attribute::read<std::int8_t>(f, "int_s"),
        Catch::Matchers::ContainsSubstring("Unable to convert field \"12345\"") &&
            Catch::Matchers::ContainsSubstring("is outside the range of representable numbers"));
  }

  SECTION("str to float lossy") {
    Attribute::write(f, "float_s",
                     std::string{"2.333333333333333481363069950020872056484222412109375"});
    CHECK(Attribute::read<float>(f, "float_s") == Catch::Approx(2.333333));
  }

  SECTION("double to int lossless") {
    CHECK(Attribute::read<std::int8_t>(f, "double") == std::int8_t(dbl));
  }

  SECTION("double to int lossy") {
    Attribute::write(f, "double", 1.1, true);
    CHECK_THROWS_WITH(Attribute::read<std::int8_t>(f, "double"),
                      Catch::Matchers::ContainsSubstring("Unable to represent value 1.1") &&
                          Catch::Matchers::ContainsSubstring("without information loss"));
  }
}

}  // namespace hictk::cooler::test::attribute
