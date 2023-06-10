// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/test/self_deleting_folder.hpp"

namespace hictk::test {
inline const SelfDeletingFolder testdir{true};                   // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::utils {
inline const auto& testdir = hictk::test::testdir;
inline const auto& datadir = hictk::test::datadir;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("utils: equal", "[equal][utils][short]") {
  const auto path1 = datadir / "cooler_test_file.cool";
  const auto path2 = datadir / "multires_cooler_test_file.mcool::/resolutions/6400000";

  const auto path3 = testdir() / "cooler_equal_test.cool";

  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path3);

  SECTION("equal") {
    CHECK(cooler::utils::equal(path1.string(), path1.string()));
    CHECK(cooler::utils::equal(path1.string(), path3.string()));
  }

  SECTION("not equal") { CHECK_FALSE(cooler::utils::equal(path1.string(), path2.string())); }
}
}  // namespace hictk::cooler::test::utils
