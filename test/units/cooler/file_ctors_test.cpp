// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: init files", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  SECTION(".cool") {
    const auto path = testdir() / "test_init.cool";
    constexpr std::uint32_t bin_size = 1000;
    std::ignore = File::create_new_cooler(path.string(), chroms, bin_size, true);
    CHECK(utils::is_cooler(path.string()));  // NOLINTNEXTLINE
    CHECK(File::open_read_only(path.string()).attributes().generated_by->find("hictk") == 0);
  }

  SECTION(".mcool") {
    const auto path = testdir() / "test_init.mcool";
    constexpr std::array<std::uint32_t, 5> resolutions{10, 20, 30, 40, 50};
    init_mcool(path.string(), resolutions.begin(), resolutions.end(), true);

    for (const auto res : resolutions) {
      std::ignore = File::create_new_cooler(
          fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), res), chroms, res);
    }

    CHECK(utils::is_multires_file(path.string()));
  }
  /*
    SECTION(".scool") {
      const auto path = (testdir() / "test_init.scool").string();
      constexpr std::array<std::string_view, 5> cell_ids{"1", "2", "3", "4", "5"};
      const std::array<std::pair<std::string_view, std::uint64_t>, 3> chroms{
          std::make_pair("chr1", 10000), std::make_pair("chr2", 5000), std::make_pair("chr3",
    1000)}; constexpr std::uint32_t bin_size = 50; init_scool(path, chroms.begin(), chroms.end(),
    cell_ids.begin(), cell_ids.end(), true);

      for (const auto id : cell_ids) {
        std::ignore =
            File::create_new_cooler(fmt::format(FMT_STRING("{}::/cells/{}"), path, id), bin_size);
      }

      CHECK(utils::is_scool_file(path));
    }
    */
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: file ctors", "[cooler][short]") {
  SECTION("default") { CHECK_NOTHROW(File{}); }

  SECTION("move #1") {
    const auto path = datadir / "cooler_test_file.cool";

    File f{};
    CHECK(!f);
    f = File::open_read_only(path.string());

    CHECK(f.chromosomes().size() == 20);
    CHECK(f.bins().size() == 26'398);
    CHECK(f.has_pixel_of_type<std::int32_t>());
  }
  SECTION("move #2") {
    const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};
    const auto path = testdir() / "move_ctor.cool";

    constexpr std::uint32_t bin_size = 1000;
    using PixelT = Pixel<std::int32_t>;
    {
      File f{};
      CHECK(!f);

      std::vector<PixelT> pixels{};
      f = File::create_new_cooler(path.string(), chroms, bin_size, true);
      const auto chr1_bins = f.bins().subset("chr1");
      for (std::uint64_t bin1_id = 0; bin1_id < chr1_bins.size(); ++bin1_id) {
        for (std::uint64_t bin2_id = bin1_id; bin2_id < chr1_bins.size(); ++bin2_id) {
          pixels.emplace_back(f.bins(), bin1_id, bin2_id,
                              static_cast<std::int32_t>(pixels.size() + 1));
        }
      }
      f.append_pixels(pixels.begin(), pixels.end(), true);
    }
  }
  SECTION("open .cool") {
    const auto path = datadir / "cooler_test_file.cool";
    const auto f = File::open_read_only(path.string());

    CHECK(f.path() == path);
    CHECK(f.uri() == path);
    CHECK(f.bin_size() == 100'000);
    CHECK(f.chromosomes().size() == 20);
    CHECK(f.bins().size() == 26'398);
    CHECK(f.has_pixel_of_type<std::int32_t>());
  }

  SECTION("open .mcool") {
    const auto path = datadir / "single_cell_cooler_test_file.scool";

    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File::open_read_only(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};
      const auto f = File::open_read_only(path.string() + suffix);

      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open .scool") {
    const auto path = datadir / "multires_cooler_test_file.mcool";
    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File::open_read_only(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[chroms, bins, pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix = "::/resolutions/400000";

      const auto f = File::open_read_only(path.string() + suffix);
      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open empty .h5") {
    const auto path = datadir / "empty_test_file.h5";
    CHECK_THROWS_WITH(File::open_read_only(path.string()),
                      Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file"));
  }

  SECTION("non existent") {
    const auto path = datadir / "cooler_test_file.cool.nonexistent";
    CHECK_THROWS_WITH(File::open_read_only(path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
  }

  SECTION("open corrupted .cool") {
    SECTION("corrupted bin table") {
      const auto path = datadir / "invalid_coolers/corrupted_bins.cool";
      CHECK_THROWS_WITH(File::open_read_only(path.string()),
                        Catch::Matchers::ContainsSubstring("Datasets have inconsistent sizes") &&
                            Catch::Matchers::ContainsSubstring("bins/chrom") &&
                            Catch::Matchers::ContainsSubstring("bins/start") &&
                            Catch::Matchers::ContainsSubstring("bins/end"));
    }

    SECTION("corrupted chrom table") {
      const auto path = datadir / "invalid_coolers/corrupted_chroms.cool";
      CHECK_THROWS_WITH(File::open_read_only(path.string()),
                        Catch::Matchers::ContainsSubstring("/chroms/name and") &&
                            Catch::Matchers::ContainsSubstring("/chroms/length shape mismatch"));
    }
  }

  SECTION("open .cool custom aprops") {
    const auto path = datadir / "cooler_test_file.cool";
    SECTION("read-once") {
      const auto f = File::open_read_only_read_once(path.string());
      CHECK(std::distance(f.begin<std::int32_t>(), f.end<std::int32_t>()) == 107041);
    }

    SECTION("read-random") {
      const auto f = File::open_read_only_random_access(path.string());
      CHECK(std::distance(f.begin<std::int32_t>(), f.end<std::int32_t>()) == 107041);
    }
  }
}

}  // namespace hictk::cooler::test::cooler_file
