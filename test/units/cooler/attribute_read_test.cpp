// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <highfive/H5File.hpp>
#include <string>

#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::attribute {

static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: attribute read", "[cooler][short]") {
  const auto path = datadir / "cooler" / "hdf5" / "test_read_attrs.h5";

  auto f = HighFive::File(path.string(), HighFive::File::ReadOnly);
  REQUIRE(f.exist("grp"));
  REQUIRE(f.exist("dst"));

  auto g = f.getGroup("grp");
  auto d = f.getDataSet("dst");

  SECTION("std::string") {
    const std::string buff{"abc"};
    SECTION("File") { CHECK(Attribute::read<std::string>(f, "std::string") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::string>(g, "std::string") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::string>(d, "std::string") == buff); }
  }

  SECTION("std::uint64_t") {
    const std::uint64_t buff{1234567890ULL};
    SECTION("File") { CHECK(Attribute::read<std::uint64_t>(f, "std::uint64_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::uint64_t>(g, "std::uint64_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::uint64_t>(d, "std::uint64_t") == buff); }
  }

  SECTION("std::int64_t") {
    const std::int64_t buff{1234567890LL};
    SECTION("File") { CHECK(Attribute::read<std::int64_t>(f, "std::int64_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::int64_t>(g, "std::int64_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::int64_t>(d, "std::int64_t") == buff); }
  }

  SECTION("double") {
    const double buff{0.123456789};
    SECTION("File") {
      CHECK_THAT(Attribute::read<double>(f, "double"), Catch::Matchers::WithinRel(buff));
    }
    SECTION("Group") {
      CHECK_THAT(Attribute::read<double>(g, "double"), Catch::Matchers::WithinRel(buff));
    }
    SECTION("Dataset") {
      CHECK_THAT(Attribute::read<double>(d, "double"), Catch::Matchers::WithinRel(buff));
    }
  }

  SECTION("std::uint32_t") {
    const std::uint32_t buff{1234567890UL};
    SECTION("File") { CHECK(Attribute::read<std::uint32_t>(f, "std::uint32_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::uint32_t>(g, "std::uint32_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::uint32_t>(d, "std::uint32_t") == buff); }
  }

  SECTION("std::int32_t") {
    const std::int32_t buff{1234567890L};
    SECTION("File") { CHECK(Attribute::read<std::int32_t>(f, "std::int32_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::int32_t>(g, "std::int32_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::int32_t>(d, "std::int32_t") == buff); }
  }

  SECTION("std::uint16_t") {
    const std::uint16_t buff{12345U};
    SECTION("File") { CHECK(Attribute::read<std::uint16_t>(f, "std::uint16_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::uint16_t>(g, "std::uint16_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::uint16_t>(d, "std::uint16_t") == buff); }
  }

  SECTION("std::int16_t") {
    const std::int16_t buff{12345};
    SECTION("File") { CHECK(Attribute::read<std::int16_t>(f, "std::int16_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::int16_t>(g, "std::int16_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::int16_t>(d, "std::int16_t") == buff); }
  }

  SECTION("std::uint8_t") {
    const std::uint8_t buff{123U};
    SECTION("File") { CHECK(Attribute::read<std::uint8_t>(f, "std::uint8_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::uint8_t>(g, "std::uint8_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::uint8_t>(d, "std::uint8_t") == buff); }
  }

  SECTION("std::int8_t") {
    const std::int8_t buff{123};
    SECTION("File") { CHECK(Attribute::read<std::int8_t>(f, "std::int8_t") == buff); }
    SECTION("Group") { CHECK(Attribute::read<std::int8_t>(g, "std::int8_t") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<std::int8_t>(d, "std::int8_t") == buff); }
  }

  SECTION("bool") {
    const bool buff{};
    SECTION("File") { CHECK(Attribute::read<bool>(f, "bool") == buff); }
    SECTION("Group") { CHECK(Attribute::read<bool>(g, "bool") == buff); }
    SECTION("Dataset") { CHECK(Attribute::read<bool>(d, "bool") == buff); }
  }

  SECTION("float") {
    const float buff{0.123456789F};
    SECTION("File") {
      CHECK_THAT(Attribute::read<float>(f, "float"), Catch::Matchers::WithinRel(buff));
    }
    SECTION("Group") {
      CHECK_THAT(Attribute::read<float>(g, "float"), Catch::Matchers::WithinRel(buff));
    }
    SECTION("Dataset") {
      CHECK_THAT(Attribute::read<float>(d, "float"), Catch::Matchers::WithinRel(buff));
    }
  }

  SECTION("std::vector") {
    SECTION("File") {
      const auto buff = Attribute::read_vector<std::int64_t>(f, "std::vector");
      REQUIRE(buff.size() == 5);
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == conditional_static_cast<std::int64_t>(i + 1));
      }
    }
    SECTION("Group") {
      const auto buff = Attribute::read_vector<std::int64_t>(g, "std::vector");
      REQUIRE(buff.size() == 5);
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == conditional_static_cast<std::int64_t>(i + 1));
      }
    }
    SECTION("Dataset") {
      const auto buff = Attribute::read_vector<std::int64_t>(d, "std::vector");
      REQUIRE(buff.size() == 5);
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == conditional_static_cast<std::int64_t>(i + 1));
      }
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::attribute
