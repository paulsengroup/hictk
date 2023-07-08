// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: version", "[cooler][short]") {
  // clang-format off
  constexpr std::array<std::uint_fast8_t, 3> ver{config::version::major(),
                                                 config::version::minor(),
                                                 config::version::patch()};
  // clang-format on

  if (config::version::suffix().empty()) {
    CHECK(HICTK_VERSION_STRING == fmt::format(FMT_STRING("{}"), fmt::join(ver, ".")));

  } else {
    CHECK(HICTK_VERSION_STRING ==
          fmt::format(FMT_STRING("{}-{}"), fmt::join(ver, "."), config::version::suffix()));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: accessors", "[cooler][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const auto f = File::open_read_only(path.string());

  SECTION("group") {
    CHECK(f.group("bins").group.getPath() == "/bins");
    CHECK_THROWS(f.group("foo"));
  }

  SECTION("dataset") {
    CHECK(f.dataset("bins/chrom").hdf5_path() == "/bins/chrom");
    CHECK_THROWS(f.dataset("/foo"));
  }

  SECTION("pixel type") {
    const auto v = f.pixel_variant();
    using T = std::int32_t;
    CHECK(std::holds_alternative<T>(v));
    CHECK(f.has_pixel_of_type<T>());

    CHECK(f.has_signed_pixels());
    CHECK_FALSE(f.has_unsigned_pixels());

    CHECK(f.has_integral_pixels());
    CHECK_FALSE(f.has_float_pixels());
  }
}

}  // namespace hictk::cooler::test::cooler_file
