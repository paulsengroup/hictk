// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tmpdir.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/test/testdir.hpp"

namespace hictk::test::common {

// NOLINTBEGIN(readability-function-cognitive-complexity)

TEST_CASE("Common: TmpDir", "[common][short]") {
  using TmpDir = internal::TmpDir;

  SECTION("Ctor") {
    SECTION("default") {
      std::filesystem::path p{};
      {
        const TmpDir tmpdir{};
        CHECK(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        p = tmpdir();
      }
      CHECK(!std::filesystem::exists(p));
    }

    SECTION("delete on destruction") {
      std::filesystem::path p{};
      {
        const TmpDir tmpdir{true};
        CHECK(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        p = tmpdir();
      }
      CHECK(!std::filesystem::exists(p));
    }

    SECTION("no delete on destruction") {
      std::filesystem::path p{};
      {
        const TmpDir tmpdir{false};
        CHECK_FALSE(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        p = tmpdir();
      }
      CHECK(std::filesystem::is_directory(p));
    }

    SECTION("from path") {
      const std::filesystem::path expected{testdir() / "unit-test-from-path"};
      {
        const TmpDir tmpdir{expected};
        CHECK(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        CHECK(tmpdir() == expected);
      }
      CHECK(!std::filesystem::exists(expected));
    }

    SECTION("from path (folder already exists)") {
      const std::filesystem::path prefix{testdir()};
      REQUIRE(std::filesystem::exists(prefix));
      CHECK_THROWS_WITH(
          TmpDir(prefix),
          Catch::Matchers::Matches("unable to use path .* as TmpDir: folder already exists"));
    }

    SECTION("from prefix (delete on destruction)") {
      const std::filesystem::path prefix{testdir() / "unit-test-from-prefix-1"};
      std::filesystem::path p{};
      REQUIRE(!std::filesystem::exists(prefix));
      std::filesystem::create_directory(prefix);
      {
        const TmpDir tmpdir{prefix, true};
        CHECK(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        CHECK(tmpdir().parent_path() == prefix);
        p = tmpdir();
      }
      CHECK(!std::filesystem::exists(p));
    }

    SECTION("from prefix (no delete on destruction)") {
      const std::filesystem::path prefix{testdir() / "unit-test-from-prefix-2"};
      std::filesystem::path p{};
      REQUIRE(!std::filesystem::exists(prefix));
      std::filesystem::create_directory(prefix);
      {
        const TmpDir tmpdir{prefix, false};
        CHECK_FALSE(tmpdir.get_delete_on_destruction());
        CHECK(!tmpdir().empty());
        CHECK(std::filesystem::is_directory(tmpdir()));
        CHECK(tmpdir().parent_path() == prefix);
        p = tmpdir();
      }
      CHECK(std::filesystem::is_directory(p));
    }

    SECTION("from prefix (prefix does not exist)") {
      const std::filesystem::path prefix{testdir() / "unit-test-from-prefix-3"};
      REQUIRE(!std::filesystem::exists(prefix));
      CHECK_THROWS_WITH(
          TmpDir(prefix, true),
          Catch::Matchers::Matches("unable to use path .* as TmpDir: path does not exists"));
    }
  }
}

// NOLINTEND(readability-function-cognitive-complexity)

}  // namespace hictk::test::common
