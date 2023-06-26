// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/dataset.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>
#include <set>

#include "hictk/cooler/group.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};                     // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::dataset {
const auto& testdir = hictk::test::testdir;
const auto& datadir = hictk::test::datadir;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset read", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      constexpr std::array<std::string_view, 3> expected{"1", "2", "3"};
      std::vector<std::string> buff{};

      Dataset{grp, "chroms/name"}.read(buff, expected.size());
      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/name"};
      std::string buff;
      dset.read(buff, 9);
      CHECK(buff == "10");
      CHECK(dset.read_last<std::string>() == "X");
      CHECK(std::get<std::string>(dset.read_last()) == "X");
    }
  }

  SECTION("numeric") {
    using T = std::int32_t;
    constexpr std::array<T, 10> expected{0,       100'000, 200'000, 300'000, 400'000,
                                         500'000, 600'000, 700'000, 800'000, 900'000};

    constexpr std::size_t nnz_expected = 107'041;
    constexpr std::int32_t sum_expected = 395'465;

    SECTION("vector<T>") {
      std::vector<T> buff{};
      std::ignore = Dataset{grp, "bins/start"}.read(buff, expected.size());

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(buff);
      CHECK(buff.size() == nnz_expected);
      CHECK(std::accumulate(buff.begin(), buff.end(), 0) == sum_expected);
    }

    SECTION("variant buff") {
      hictk::internal::VariantBuffer vbuff{std::size_t(0), 0.0};
      std::ignore = Dataset{grp, "bins/start"}.read(vbuff, expected.size());
      const auto& buff = vbuff.get<T>();
      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(vbuff);
      CHECK(vbuff.size<T>() == nnz_expected);
      CHECK(std::accumulate(vbuff.begin<T>(), vbuff.end<T>(), 0) == sum_expected);
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/length"};
      std::uint64_t buff{};
      dset.read(buff, 2);
      CHECK(buff == 159'599'783);

      CHECK(dset.read_last<std::int32_t>() == 166'650'296);
      CHECK(std::get<std::int32_t>(dset.read_last()) == 166'650'296);
    }

    SECTION("enum") { CHECK(Dataset{grp, "bins/chrom"}.read<std::uint32_t>(0) == 0); }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset accessors", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  Dataset dset{grp, "chroms/name"};

  CHECK(dset.size() == 20);

  CHECK(dset.get().getFile().getName() == path.string());

  CHECK(dset.file_name() == path.string());
  CHECK(dset.uri() == fmt::format(FMT_STRING("{}::/chroms/name"), path.string()));
  CHECK(dset.hdf5_path() == "/chroms/name");
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset linear iteration", "[dataset][long]") {
  const auto path = datadir / "cooler_test_file.cool";

  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
  const Dataset dset(grp, "/pixels/count");

  std::vector<std::uint32_t> pixel_buff;
  dset.read_all(pixel_buff);
  REQUIRE(pixel_buff.size() == 107'041);

  SECTION("forward") {
    auto it = dset.begin<std::uint32_t>();
    auto last_pixel = dset.end<std::uint32_t>();
    REQUIRE(std::distance(it, last_pixel) == 107'041);

    for (const auto& expected : pixel_buff) {
      REQUIRE(it != last_pixel);
      CHECK(*it++ == expected);
    }
    CHECK(it == last_pixel);
  }

  SECTION("backward") {
    auto it = dset.end<std::uint32_t>();
    auto first_pixel = dset.begin<std::uint32_t>();

    REQUIRE(std::distance(first_pixel, it) == 107'041);

    for (std::size_t i = pixel_buff.size(); i != 0; --i) {
      REQUIRE(it != first_pixel);
      CHECK(*(--it) == pixel_buff[i - 1]);
    }

    CHECK(it == first_pixel);
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset random iteration", "[dataset][medium]") {
  const auto path = testdir() / "dataset_iterator_random.h5";

  const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
  Dataset dset(grp, "int", std::uint64_t{});

  std::random_device rd;
  std::mt19937_64 rand_eng{rd()};

  constexpr std::size_t N = 5'000'000;
  std::vector<std::uint64_t> buff(N);
  std::generate(buff.begin(), buff.end(), [&]() { return rand_eng(); });
  dset.append(buff);
  REQUIRE(dset.size() == N);

  SECTION("operator -/+") {
    auto first = dset.begin<std::uint64_t>();
    auto last = dset.end<std::uint64_t>();
    for (std::size_t i = 0; i < 100; ++i) {
      const auto j = std::uniform_int_distribution<std::uint64_t>{0, N - 1}(rand_eng);

      CHECK(*(first + j) == buff[j]);
      CHECK(*(last - j) == buff[N - j]);
    }
  }

  SECTION("subsequent calls to operator+=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.begin<std::uint64_t>();
      auto last = dset.end<std::uint64_t>();
      std::size_t j = 0;

      while (first < last) {
        CHECK(*first == buff[j]);

        const auto step = std::uniform_int_distribution<std::size_t>{
            0, (std::min)(std::size_t(500), buff.size() - j)}(rand_eng);
        first += step;
        j += step;
      }
    }
  }

  SECTION("subsequent calls to operator-=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.end<std::uint64_t>() - 1;
      auto last = dset.begin<std::uint64_t>();
      std::size_t j = buff.size() - 1;

      while (first > last) {
        CHECK(*first == buff[j]);

        const auto step = std::uniform_int_distribution<std::size_t>{
            0, (std::min)(std::size_t(500), j)}(rand_eng);

        first -= step;
        j -= step;
      }
    }
  }
}

TEST_CASE("Cooler: dataset large read/write", "[dataset][long]") {
  const auto path = testdir() / "test_dataset_large_rw.h5";

  constexpr std::uint64_t seed{4195331987557451569};
  constexpr std::size_t N = 5'000'000;
  {
    const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    std::vector<std::uint8_t> buff(1'000'000);
    std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
    while (dset.size() != N) {
      std::generate(buff.begin(), buff.end(),
                    [&]() { return static_cast<std::uint8_t>(rand_eng()); });
      dset.append(buff);
    }
    CHECK(dset.size() == N);
  }

  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
  const Dataset dset{grp, "int"};
  REQUIRE(dset.size() == N);

  std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
  std::for_each(dset.begin<std::uint8_t>(), dset.end<std::uint8_t>(),
                [&](const auto& n) { CHECK(n == static_cast<std::uint8_t>(rand_eng())); });
}

TEST_CASE("Cooler: dataset attributes", "[dataset][short]") {
  SECTION("read") {
    const auto path = datadir / "test_read_attrs.h5";

    const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
    const Dataset dset{grp, "dst"};

    CHECK(dset.has_attribute("std::string"));
    CHECK(dset.read_attribute<std::string>("std::string") == "abc");
    auto attr = dset.read_attribute("std::string");
    CHECK(std::holds_alternative<std::string>(attr));
    CHECK_THROWS(dset.read_attribute("invalid"));
    attr = dset.read_attribute("invalid", true);
    CHECK(std::holds_alternative<std::monostate>(attr));
  }

  SECTION("write") {
    const auto path = testdir() / "test_dataset_write_attr.h5";

    const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    dset.write_attribute("attr", 123);
    CHECK(dset.read_attribute<std::int32_t>("attr") == 123);

    CHECK_THROWS(dset.write_attribute("attr", -1, false));
    dset.write_attribute("attr", -1, true);
    CHECK(dset.read_attribute<std::int32_t>("attr") == -1);
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

}  // namespace hictk::cooler::test::dataset
