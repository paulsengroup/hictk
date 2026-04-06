// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>

// clang-format off
#include "hictk/git.hpp"
#include "hictk/license.hpp"
#include "hictk/version.hpp"
// clang-format on

namespace hictk::test::common {

// NOLINTBEGIN(readability-function-cognitive-complexity)

TEST_CASE("Common: git metadata", "[common][short]") {
  namespace git = config::git;
  if constexpr (git::state_available()) {
    SECTION("state available") {
      CHECK(git::head_sha1().size() == 40);
      CHECK_FALSE(git::author_name().empty());
      CHECK_FALSE(git::author_email().empty());
      CHECK_FALSE(git::commit_date().empty());
      CHECK_FALSE(git::commit_subject().empty());
      CHECK_FALSE(git::describe().empty());
      CHECK_FALSE(git::branch().empty());
    }
  } else {
    SECTION("state not available") {
      CHECK(git::head_sha1() == "unknown");
      CHECK_FALSE(git::is_dirty());
      CHECK(git::author_name() == "unknown");
      CHECK(git::author_email() == "unknown");
      CHECK(git::commit_date() == "unknown");
      CHECK(git::commit_subject() == "unknown");
      CHECK(git::commit_body() == "unknown");
      CHECK(git::describe() == "unknown");
      CHECK(git::branch() == "unknown");
      CHECK(git::tag() == "unknown");
    }
  }
}

TEST_CASE("Common: license metadata", "[common][short]") {
  const std::string license{config::license::license};
  REQUIRE_FALSE(license.empty());
  CHECK_THAT(license, Catch::Matchers::StartsWith("MIT License"));
  CHECK_THAT(license, Catch::Matchers::ContainsSubstring("Copyright"));
}

TEST_CASE("Common: version metadata", "[common][short]") {
  namespace version = config::version;
  const auto suffix_avail = !version::suffix().empty();

  SECTION("string") {
    auto version_str{
        fmt::format(FMT_STRING("{}.{}.{}"), version::major(), version::minor(), version::patch())};

    if (suffix_avail) {
      version_str += fmt::format(FMT_STRING("-{}"), version::suffix());
      CHECK_THAT(std::string{version::str()}, Catch::Matchers::StartsWith(version_str));
    } else {
      CHECK_THAT(std::string{version::str()}, Catch::Matchers::Equals(version_str));
    }
  }

  SECTION("string long") {
    auto version_str{fmt::format(FMT_STRING("hictk-v{}.{}.{}"), version::major(), version::minor(),
                                 version::patch())};
    if (suffix_avail) {
      version_str += fmt::format(FMT_STRING("-{}"), version::suffix());
      CHECK_THAT(std::string{version::str_long()}, Catch::Matchers::StartsWith(version_str));
    } else {
      CHECK_THAT(std::string{version::str_long()}, Catch::Matchers::Equals(version_str));
    }
  }

  SECTION("global constants") {
    CHECK_THAT(std::string{HICTK_VERSION_STRING},
               Catch::Matchers::Equals(std::string{version::str()}));
    CHECK_THAT(std::string{HICTK_VERSION_STRING_LONG},
               Catch::Matchers::Equals(std::string{version::str_long()}));
  }
}

// NOLINTEND(readability-function-cognitive-complexity)

}  // namespace hictk::test::common
