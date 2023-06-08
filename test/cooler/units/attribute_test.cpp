// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/attribute.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <highfive/H5File.hpp>
#include <string>

#include "coolerpp/internal/suppress_warnings.hpp"
#include "coolerpp/test/self_deleting_folder.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};                   // NOLINT(cert-err58-cpp)
static inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::attribute {

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Attribute: write", "[cooler][short]") {
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Attribute: read", "[cooler][short]") {
  const auto path = datadir / "test_read_attrs.h5";

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Attribute: read - test numeric conversions", "[cooler][short]") {
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

}  // namespace coolerpp::test::attribute
