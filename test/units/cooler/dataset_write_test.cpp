// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/test/testdir.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk::cooler::test::dataset {

static const auto& testdir = hictk::test::testdir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset write", "[dataset][short]") {
  const auto path = testdir() / "test_dataset_write.cool";
  const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      using BuffT = std::vector<std::string>;
      const BuffT expected{"s1", "this_is_a_relatively_long_string"};

      Dataset{grp, "str", expected.back()}.write(expected, 0, true);
      const auto buff = Dataset{grp, "str"}.read_all<BuffT>();

      REQUIRE(buff.size() == 2);

      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("iterator") {
      const std::set<std::string> expected{{"a", "b", "c"}};
      Dataset{grp, "str", "a"}.write(expected.begin(), expected.end(), 0, true);
      const auto buff = Dataset{grp, "str"}.read_all<std::vector<std::string>>();

      REQUIRE(buff.size() == 3);

      for (const auto& i : buff) {
        CHECK(expected.count(i) == 1);
      }
    }

    SECTION("pointer") {
      const std::array<std::string, 3> expected{"a", "b", "c"};
      Dataset{grp, "str", "a"}.write(expected.begin(), expected.end(), 0, true);
      const auto buff = Dataset{grp, "str"}.read_all<std::vector<std::string>>();

      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      const std::string buff = "test";
      Dataset{grp, "str", buff}.write(buff, 3, true);

      CHECK(Dataset{grp, "str"}.read<std::string>(0).empty());
      CHECK(Dataset{grp, "str"}.read<std::string>(3) == buff);
    }
  }

  SECTION("numeric") {
    using T = double;
    using BuffT = std::vector<T>;
    const BuffT expected{0.1, 0.2, 0.3};

    SECTION("vector<N>") {
      Dataset{grp, "num", T{}}.write(expected, 0, true);
      const auto buff = Dataset{grp, "num"}.read_all<BuffT>();

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("variant buff") {
      const hictk::internal::VariantBuffer vexpected{expected};
      Dataset{grp, "num", T{}}.write(vexpected, 0, true);

      const auto vbuff = Dataset{grp, "num"}.read_all();
      const auto& buff = vbuff.get<T>();

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      Dataset{grp, "num", T{}}.write(7.0, 5, true);
      REQUIRE(Dataset{grp, "num"}.size() == 6);

      CHECK(Dataset{grp, "num"}.read<T>(0) == 0.0);
      CHECK(Dataset{grp, "num"}.read<T>(5) == 7.0);

      const auto vbuff = Dataset{grp, "num"}.read(5);
      CHECK(std::get<T>(vbuff) == 7.0);
    }
  }

  SECTION("out of bound access") {
    CHECK_THROWS_WITH((Dataset{grp, "num", int{}}.write(1, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("which is empty"));

    Dataset{grp, "num"}.resize(10);
    CHECK_THROWS_WITH((Dataset{grp, "num"}.write(1, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("past the end"));

    CHECK_THROWS_WITH((Dataset{grp, "num"}.write(std::vector<int>{1, 2, 3}, 100U, false)),
                      Catch::Matchers::ContainsSubstring("attempt to access") &&
                          Catch::Matchers::ContainsSubstring("past the end"));
  }
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::dataset
