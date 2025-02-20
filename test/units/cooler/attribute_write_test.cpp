// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <highfive/H5File.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/attribute.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::attribute {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
template <typename H5Obj, typename T>
static void compare_attribute(H5Obj& obj, std::string_view key, const T& expected) {
  static_assert(std::is_same_v<T, std::string> || std::is_fundamental_v<T>);
  T buff{};
  obj.getAttribute(std::string{key}).read(buff);

  CHECK(expected == buff);
}

template <typename H5Obj, typename T>
static void compare_attribute(H5Obj& obj, std::string_view key, const std::vector<T>& expected) {
  std::vector<T> buff{};
  obj.getAttribute(std::string{key}).read(buff);
  REQUIRE(expected.size() == buff.size());
  for (std::size_t i = 0; i < buff.size(); ++i) {
    CHECK(expected[i] == buff[i]);
  }
}

TEST_CASE("Cooler: attribute write", "[cooler][short]") {
  const auto path = testdir() / "test_write_attrs.cool";

  auto f = HighFive::File(path.string(), HighFive::File::Truncate);
  auto g = f.createGroup("grp");
  auto d = f.createDataSet("dst", std::string{});

  SECTION("std::string") {
    const std::string buff{"abc"};
    SECTION("File") {
      Attribute::write(f, "std::string", buff);
      compare_attribute(f, "std::string", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::string", buff);
      compare_attribute(g, "std::string", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::string", buff);
      compare_attribute(d, "std::string", buff);
    }
  }

  SECTION("std::uint64_t") {
    const std::uint64_t buff{1234567890ULL};
    SECTION("File") {
      Attribute::write(f, "std::uint64_t", buff);
      compare_attribute(f, "std::uint64_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::uint64_t", buff);
      compare_attribute(g, "std::uint64_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::uint64_t", buff);
      compare_attribute(d, "std::uint64_t", buff);
    }
  }

  SECTION("std::int64_t") {
    const std::int64_t buff{1234567890LL};
    SECTION("File") {
      Attribute::write(f, "std::int64_t", buff);
      compare_attribute(f, "std::int64_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::int64_t", buff);
      compare_attribute(g, "std::int64_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::int64_t", buff);
      compare_attribute(d, "std::int64_t", buff);
    }
  }

  SECTION("double") {
    const double buff{0.123456789};
    SECTION("File") {
      Attribute::write(f, "double", buff);
      compare_attribute(f, "double", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "double", buff);
      compare_attribute(g, "double", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "double", buff);
      compare_attribute(d, "double", buff);
    }
  }

  SECTION("std::uint32_t") {
    const std::uint32_t buff{1234567890UL};
    SECTION("File") {
      Attribute::write(f, "std::uint32_t", buff);
      compare_attribute(f, "std::uint32_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::uint32_t", buff);
      compare_attribute(g, "std::uint32_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::uint32_t", buff);
      compare_attribute(d, "std::uint32_t", buff);
    }
  }

  SECTION("std::int32_t") {
    const std::int32_t buff{1234567890L};
    SECTION("File") {
      Attribute::write(f, "std::int32_t", buff);
      compare_attribute(f, "std::int32_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::int32_t", buff);
      compare_attribute(g, "std::int32_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::int32_t", buff);
      compare_attribute(d, "std::int32_t", buff);
    }
  }

  SECTION("std::uint16_t") {
    const std::uint16_t buff{12345U};
    SECTION("File") {
      Attribute::write(f, "std::uint16_t", buff);
      compare_attribute(f, "std::uint16_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::uint16_t", buff);
      compare_attribute(g, "std::uint16_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::uint16_t", buff);
      compare_attribute(d, "std::uint16_t", buff);
    }
  }

  SECTION("std::int16_t") {
    const std::int16_t buff{12345};
    SECTION("File") {
      Attribute::write(f, "std::int16_t", buff);
      compare_attribute(f, "std::int16_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::int16_t", buff);
      compare_attribute(g, "std::int16_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::int16_t", buff);
      compare_attribute(d, "std::int16_t", buff);
    }
  }

  SECTION("std::uint8_t") {
    const std::uint8_t buff{123U};
    SECTION("File") {
      Attribute::write(f, "std::uint8_t", buff);
      compare_attribute(f, "std::uint8_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::uint8_t", buff);
      compare_attribute(g, "std::uint8_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::uint8_t", buff);
      compare_attribute(d, "std::uint8_t", buff);
    }
  }

  SECTION("std::int8_t") {
    const std::int8_t buff{123};
    SECTION("File") {
      Attribute::write(f, "std::int8_t", buff);
      compare_attribute(f, "std::int8_t", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::int8_t", buff);
      compare_attribute(g, "std::int8_t", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::int8_t", buff);
      compare_attribute(d, "std::int8_t", buff);
    }
  }

  SECTION("bool") {
    const bool buff{};
    SECTION("File") {
      Attribute::write(f, "bool", buff);
      compare_attribute(f, "bool", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "bool", buff);
      compare_attribute(g, "bool", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "bool", buff);
      compare_attribute(d, "bool", buff);
    }
  }

  SECTION("long double") {
    const long double buff{0.123456789L};
    SECTION("File") {
      Attribute::write(f, "long double", buff);
      compare_attribute(f, "long double", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "long double", buff);
      compare_attribute(g, "long double", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "long double", buff);
      compare_attribute(d, "long double", buff);
    }
  }

  SECTION("float") {
    const float buff{0.123456789F};
    SECTION("File") {
      Attribute::write(f, "float", buff);
      compare_attribute(f, "float", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "float", buff);
      compare_attribute(g, "float", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "float", buff);
      compare_attribute(d, "float", buff);
    }
  }

  SECTION("std::vector") {
    const std::vector<std::int32_t> buff{1, 2, 3};
    SECTION("File") {
      Attribute::write(f, "std::vector", buff);
      compare_attribute(f, "std::vector", buff);
    }

    SECTION("Group") {
      Attribute::write(g, "std::vector", buff);
      compare_attribute(g, "std::vector", buff);
    }

    SECTION("Dataset") {
      Attribute::write(d, "std::vector", buff);
      compare_attribute(d, "std::vector", buff);
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::attribute
