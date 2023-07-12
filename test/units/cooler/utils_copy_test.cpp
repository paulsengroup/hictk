// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};                     // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::utils {
inline const auto& testdir = hictk::test::testdir;
inline const auto& datadir = hictk::test::datadir;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils copy", "[merge][utils][long]") {
  SECTION("cooler -> cooler") {
    const auto src = datadir / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_copy_001.cool";

    cooler::utils::copy(src.string(), dest.string());
    CHECK(cooler::utils::equal(src.string(), dest.string()));
  }

  SECTION("cooler -> mcool") {
    const auto src = datadir / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_copy_002.mcool";
    const auto dest_uri = fmt::format(FMT_STRING("{}::/resolutions/1000"), dest.string());

    std::vector<std::uint32_t> resolutions{1000};
    cooler::init_mcool(dest.string(), resolutions.begin(), resolutions.end());

    cooler::utils::copy(src.string(), dest_uri);
    CHECK(cooler::utils::equal(src.string(), dest_uri));
  }

  SECTION("mcool -> cooler") {
    const auto src = datadir / "multires_cooler_test_file.mcool::/resolutions/100000";
    const auto dest = testdir() / "cooler_copy_003.cool";

    cooler::utils::copy(src.string(), dest.string());
    CHECK(cooler::utils::equal(src.string(), dest.string()));
  }
}

}  // namespace hictk::cooler::test::utils
