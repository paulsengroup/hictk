// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "coolerpp/test/self_deleting_folder.hpp"
#include "coolerpp/utils.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};            // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::index {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("utils: equal", "[equal][utils][short]") {
  const auto path1 = datadir / "cooler_test_file.cool";
  const auto path2 = datadir / "multires_cooler_test_file.mcool::/resolutions/6400000";

  const auto path3 = testdir() / "cooler_equal_test.cool";

  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path3);

  SECTION("equal") {
    CHECK(utils::equal(path1.string(), path1.string()));
    CHECK(utils::equal(path1.string(), path3.string()));
  }

  SECTION("not equal") { CHECK_FALSE(utils::equal(path1.string(), path2.string())); }
}
}  // namespace coolerpp::test::index
