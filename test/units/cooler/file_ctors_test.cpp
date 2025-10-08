// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <tuple>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::cooler_file {

static const auto& datadir = hictk::test::datadir;
static const auto& testdir = hictk::test::testdir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: create files", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  SECTION("fixed bins") {
    const auto path = testdir() / "test_init_fixed_bins.cool";
    constexpr std::uint32_t bin_size = 1000;
    std::ignore = File::create(path.string(), chroms, bin_size, true);
    CHECK(utils::is_cooler(path.string()));  // NOLINTNEXTLINE
    CHECK(File(path.string()).attributes().generated_by->find("hictk") == 0);
    CHECK(File(path.string()).attributes().bin_type == BinTable::Type::fixed);
  }

  SECTION("variable bins") {
    const auto path = testdir() / "test_init_variable_bins.cool";

    const Chromosome chrom1{0, "chr1", 32};
    const Chromosome chrom2{1, "chr2", 32};

    const std::vector<std::uint32_t> start_pos{0, 8, 15, 23, 0, 5, 10, 26};
    const std::vector<std::uint32_t> end_pos{8, 15, 23, 32, 5, 10, 26, 32};

    const BinTable table({chrom1, chrom2}, start_pos, end_pos);

    std::ignore = File::create(path.string(), table, true);
    CHECK(utils::is_cooler(path.string()));  // NOLINTNEXTLINE
    CHECK(File(path.string()).attributes().generated_by->find("hictk") == 0);
    CHECK(File(path.string()).attributes().bin_type == BinTable::Type::variable);
  }

  // NOLINTBEGIN(*-pointer-arithmetic)
  SECTION("append pixels") {
    SECTION("valid") {
      const auto path = testdir() / "test_init_append_pixels_valid.cool";
      constexpr std::uint32_t bin_size = 1000;
      auto clr = File::create(path.string(), chroms, bin_size, true);

      const ThinPixel<std::int32_t> p1{0, 0, 1};
      const Pixel p2{clr.bins(), ThinPixel<std::int32_t>{0, 1, 1}};

      clr.append_pixels(&p1, &p1 + 1, true);
      CHECK(clr.attributes().nnz == 1);
      clr.append_pixels(&p2, &p2 + 1, true);
      CHECK(clr.attributes().nnz == 2);
    }

    SECTION("invalid") {
      const auto path = testdir() / "test_init_append_pixels_invalid.cool";
      constexpr std::uint32_t bin_size = 1000;
      auto clr = File::create(path.string(), chroms, bin_size, true);
      const BinTable invalid_bins{clr.chromosomes(), bin_size / 2};

      SECTION("invalid count") {
        const ThinPixel<std::int32_t> p1{0, 0, 0};
        const Pixel p2{clr.bins(), ThinPixel<std::int32_t>{0, 1, 0}};

        CHECK_THROWS_WITH(clr.append_pixels(&p1, &p1 + 1, true),
                          Catch::Matchers::ContainsSubstring("found a pixel of value 0"));
        CHECK_THROWS_WITH(clr.append_pixels(&p2, &p2 + 1, true),
                          Catch::Matchers::ContainsSubstring("found a pixel of value 0"));
      }

      SECTION("invalid chrom1") {
        const Chromosome chrom{2, "chr3", 10000};
        const Bin bin1{0, 0, chrom, 0, bin_size};
        const Bin bin2{0, 0, clr.chromosomes().at("chr1"), 0, bin_size};
        const Pixel p{bin1, bin2, 1};

        CHECK_THROWS_WITH(clr.append_pixels(&p, &p + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid chromosome id"));
      }

      SECTION("invalid chrom2") {
        const Chromosome chrom{2, "chr3", 10000};
        const Bin bin1{0, 0, clr.chromosomes().at("chr1"), 0, bin_size};
        const Bin bin2{0, 0, chrom, 0, bin_size};
        const Pixel p{bin1, bin2, 1};

        CHECK_THROWS_WITH(clr.append_pixels(&p, &p + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid chromosome id"));
      }

      SECTION("invalid bin1_id") {
        const ThinPixel<std::int32_t> p1{16, 16, 1};
        const Pixel p2{invalid_bins, 16, 16, 1};

        CHECK_THROWS_WITH(clr.append_pixels(&p1, &p1 + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid bin id"));
        CHECK_THROWS_WITH(clr.append_pixels(&p2, &p2 + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid bin id"));
      }

      SECTION("invalid bin2_id") {
        const ThinPixel<std::int32_t> p1{0, 16, 1};
        const Pixel p2{invalid_bins, 0, 16, 1};

        CHECK_THROWS_WITH(clr.append_pixels(&p1, &p1 + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid bin id"));
        CHECK_THROWS_WITH(clr.append_pixels(&p2, &p2 + 1, true),
                          Catch::Matchers::ContainsSubstring("invalid bin id"));
      }

      SECTION("lower triangle") {
        const ThinPixel<std::int32_t> p1{1, 0, 1};
        const Pixel p2{clr.bins(), 1, 0, 1};

        CHECK_THROWS_WITH(clr.append_pixels(&p1, &p1 + 1, true),
                          Catch::Matchers::ContainsSubstring("bin1_id is greater than bin2_id"));
        CHECK_THROWS_WITH(clr.append_pixels(&p2, &p2 + 1, true),
                          Catch::Matchers::ContainsSubstring("bin1_id is greater than bin2_id"));
      }

      SECTION("out of order chunks") {
        const ThinPixel<std::int32_t> p1{1, 2, 1};
        const ThinPixel<std::int32_t> p2{1, 1, 1};
        const ThinPixel<std::int32_t> p3{0, 0, 1};
        const Pixel p4{clr.bins(), p2};
        const Pixel p5{clr.bins(), p3};

        clr.append_pixels(&p1, &p1 + 1, true);

        CHECK_THROWS_WITH(clr.append_pixels(&p2, &p2 + 1, true),
                          Catch::Matchers::ContainsSubstring("new pixel") &&
                              Catch::Matchers::ContainsSubstring("is located upstream"));
        CHECK_THROWS_WITH(clr.append_pixels(&p3, &p3 + 1, true),
                          Catch::Matchers::ContainsSubstring("new pixel") &&
                              Catch::Matchers::ContainsSubstring("is located upstream"));
        CHECK_THROWS_WITH(clr.append_pixels(&p4, &p4 + 1, true),
                          Catch::Matchers::ContainsSubstring("new pixel") &&
                              Catch::Matchers::ContainsSubstring("is located upstream"));
        CHECK_THROWS_WITH(clr.append_pixels(&p5, &p5 + 1, true),
                          Catch::Matchers::ContainsSubstring("new pixel") &&
                              Catch::Matchers::ContainsSubstring("is located upstream"));
      }

      SECTION("unsorted") {
        const std::array<ThinPixel<std::int32_t>, 2> p1{ThinPixel<std::int32_t>{0, 1, 1},
                                                        ThinPixel<std::int32_t>{0, 0, 1}};
        const std::array<Pixel<std::int32_t>, 2> p2{Pixel{clr.bins(), p1.front()},
                                                    Pixel{clr.bins(), p1.back()}};

        CHECK_THROWS_WITH(clr.append_pixels(p1.begin(), p1.end(), true),
                          Catch::Matchers::ContainsSubstring("pixels are not sorted"));
        CHECK_THROWS_WITH(clr.append_pixels(p2.begin(), p2.end(), true),
                          Catch::Matchers::ContainsSubstring("pixels are not sorted"));
      }
    }
  }
  // NOLINTEND(*-pointer-arithmetic)
}

TEST_CASE("Cooler: file ctors", "[cooler][short]") {
  SECTION("default") {
    CHECK_NOTHROW(File{});
    CHECK(File{}.path() == "");
    CHECK(File{}.hdf5_path() == "");
    CHECK(File{}.uri() == "");
  }

  SECTION("move #1") {
    const auto path = datadir / "cooler" / "cooler_test_file.cool";

    File f{};
    CHECK(!f);
    f = File(path.string());

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
      f = File::create(path.string(), chroms, bin_size, true);
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
  SECTION("open .cool (fixed bin size)") {
    const auto path = datadir / "cooler" / "cooler_test_file.cool";
    const File f(path.string());

    CHECK(f.path() == path);
    CHECK(f.uri() == path);
    CHECK(f.resolution() == 100'000);
    CHECK(f.chromosomes().size() == 20);
    CHECK(f.bins().size() == 26'398);
    CHECK(f.has_pixel_of_type<std::int32_t>());
  }

  SECTION("open .cool (variable bin size)") {
    const auto path = datadir / "cooler" / "cooler_variable_bins_test_file.cool";
    const File f(path.string());

    CHECK(f.path() == path);
    CHECK(f.uri() == path);
    CHECK(f.resolution() == 0);
    CHECK(f.chromosomes().size() == 2);
    CHECK(f.bins().size() == 8);
    CHECK(f.has_pixel_of_type<std::int32_t>());
  }

  SECTION("open .cool (storage-mode=square)") {
    const auto path =
        datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
    const File f(path.string());

    CHECK(f.uri() == path.string());
    CHECK(f.resolution() == 1000);
    CHECK(f.chromosomes().size() == 10);
    CHECK(f.bins().size() == 3000);
    CHECK(f.has_pixel_of_type<std::int32_t>());
    CHECK(f.attributes().storage_mode == "square");
  }

  SECTION("open .scool") {
    const auto path = datadir / "cooler" / "single_cell_cooler_test_file.scool";

    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};
      const File f(path.string() + suffix);

      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open .mcool") {
    const auto path = datadir / "cooler" / "multires_cooler_test_file.mcool";
    SECTION("missing suffix") {
      CHECK_THROWS_WITH(
          File(path.string()),
          Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file") &&
              Catch::Matchers::ContainsSubstring("missing_groups=[chroms, bins, pixels, indexes]"));
    }

    SECTION("with suffix") {
      constexpr auto suffix = "::/resolutions/400000";

      const File f(path.string() + suffix);
      CHECK(f.path() == path);
      CHECK(f.uri() == path.string() + suffix);
    }
  }

  SECTION("open empty .h5") {
    const auto path = datadir / "cooler" / "hdf5" / "empty_test_file.h5";
    CHECK_THROWS_WITH(File(path.string()),
                      Catch::Matchers::ContainsSubstring("does not look like a valid Cooler file"));
  }

  SECTION("non existent") {
    const auto path = datadir / "cooler_test_file.cool.nonexistent";
    CHECK_THROWS_WITH(File(path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
  }

  SECTION("open corrupted .cool") {
    SECTION("corrupted bin table") {
      const auto path = datadir / "cooler" / "invalid" / "corrupted_bins.cool";
      CHECK_THROWS_WITH(File(path.string()),
                        Catch::Matchers::ContainsSubstring("Datasets have inconsistent sizes") &&
                            Catch::Matchers::ContainsSubstring("bins/chrom") &&
                            Catch::Matchers::ContainsSubstring("bins/start") &&
                            Catch::Matchers::ContainsSubstring("bins/end"));
    }

    SECTION("corrupted chrom table") {
      const auto path = datadir / "cooler" / "invalid" / "corrupted_chroms.cool";
      CHECK_THROWS_WITH(File(path.string()),
                        Catch::Matchers::ContainsSubstring("/chroms/name and") &&
                            Catch::Matchers::ContainsSubstring("/chroms/length shape mismatch"));
    }
  }

  SECTION("open .cool custom aprops") {
    const auto path = datadir / "cooler" / "cooler_test_file.cool";
    SECTION("read-once") {
      const auto f = File::open_read_once(path.string());
      CHECK(std::distance(f.begin<std::int32_t>(), f.end<std::int32_t>()) == 107041);
    }

    SECTION("read-random") {
      const auto f = File::open_random_access(path.string());
      CHECK(std::distance(f.begin<std::int32_t>(), f.end<std::int32_t>()) == 107041);
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::cooler_file
