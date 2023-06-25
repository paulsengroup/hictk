// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>

#include "hictk/cooler.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};                     // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::cooler_file {

const auto& testdir = hictk::test::testdir;
const auto& datadir = hictk::test::datadir;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: version", "[cooler][short]") {
  // clang-format off
  constexpr std::array<std::uint_fast8_t, 3> ver{config::version::major(),
                                                 config::version::minor(),
                                                 config::version::patch()};
  // clang-format on

  if (config::version::suffix().empty()) {
    CHECK(HICTK_VERSION_STRING == fmt::format(FMT_STRING("{}"), fmt::join(ver, ".")));

  } else {
    CHECK(HICTK_VERSION_STRING ==
          fmt::format(FMT_STRING("{}-{}"), fmt::join(ver, "."), config::version::suffix()));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: format checking", "[cooler][short]") {
  SECTION("test .cool") {
    const auto path = datadir / "cooler_test_file.cool";
    CHECK(utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test .mcool") {
    const auto path = datadir / "multires_cooler_test_file.mcool";
    constexpr auto suffix{"::/resolutions/400000"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test .scool") {
    const auto path = datadir / "single_cell_cooler_test_file.scool";
    constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test with empty .h5 file") {
    const auto path = datadir / "empty_test_file.h5";
    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test with nonexistent file") {
    const auto invalid_path = datadir / "void.nonexistent";
    CHECK_THROWS_WITH(utils::is_cooler(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
    CHECK_THROWS_WITH(utils::is_multires_file(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
    CHECK_THROWS_WITH(utils::is_scool_file(invalid_path.string()),
                      Catch::Matchers::ContainsSubstring("Unable to open file"));
  }

  SECTION("test corrupted .cool") {
    SECTION("missing format attribute") {
      const auto path = datadir / "invalid_coolers/missing_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
    SECTION("invalid format attribute") {
      const auto path = datadir / "invalid_coolers/invalid_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
  }

  SECTION("test corrupted .mcool") {
    // This file is missing group /resolutions/400000/pixels
    const auto path = datadir / "invalid_coolers/missing_pixels_group.mcool";
    const auto status = utils::is_multires_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_multires_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());

    REQUIRE(status.invalid_resolutions.size() == 1);
    const auto& invalid_res = status.invalid_resolutions.front();

    const auto corrupted_uri_expected =
        fmt::format(FMT_STRING("{}::/resolutions/400000"), path.string());
    CHECK(invalid_res.uri == corrupted_uri_expected);
    CHECK(!invalid_res.is_cooler);
    REQUIRE(invalid_res.missing_groups.size() == 1);
    CHECK(invalid_res.missing_groups.front() == "pixels");
  }

  SECTION("test corrupted .scool") {
    // In this file, the number of groups under /cells and number of cells from ncells attribute
    // mismatch
    const auto path = datadir / "invalid_coolers/invalid_ncells_attribute.scool";
    const auto status = utils::is_scool_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_scool_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());
    CHECK(status.unexpected_number_of_cells);
    CHECK(status.invalid_cells.empty());
  }
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: accessors", "[cooler][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const auto f = File::open_read_only(path.string());

  SECTION("group") {
    CHECK(f.group("bins").group.getPath() == "/bins");
    CHECK_THROWS(f.group("foo"));
  }

  SECTION("dataset") {
    CHECK(f.dataset("bins/chrom").hdf5_path() == "/bins/chrom");
    CHECK_THROWS(f.dataset("foo"));
  }

  SECTION("pixel type") {
    const auto v = f.pixel_variant();
    using T = std::int32_t;
    CHECK(std::holds_alternative<T>(v));
    CHECK(f.has_pixel_of_type<T>());

    CHECK(f.has_signed_pixels());
    CHECK_FALSE(f.has_unsigned_pixels());

    CHECK(f.has_integral_pixels());
    CHECK_FALSE(f.has_float_pixels());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: init files", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  SECTION(".cool") {
    const auto path = testdir() / "test_init.cool";
    constexpr std::uint32_t bin_size = 1000;
    std::ignore = File::create_new_cooler(path.string(), chroms, bin_size, true);
    CHECK(utils::is_cooler(path.string()));
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
TEST_CASE("Cooler: sentinel attribute", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  const auto path = testdir() / "test_sentinel_attr.cool";
  constexpr std::uint32_t bin_size = 1000;
  auto f = File::create_new_cooler(path.string(), chroms, bin_size, true);

  SECTION("Read-only") {
    const auto path1 = datadir / "cooler_test_file.cool";
    const auto f1 = File::open_read_only(path1.string());
    CHECK(Attribute::read<std::uint8_t>(f1.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);
    f.close();
    f = File::open_read_only(path.string());
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create (file was not closed properly)") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);

    CHECK_THROWS(f = File::open_read_only(path.string()));
    CHECK_THROWS(f = File::create_new_cooler(path.string(), chroms, bin_size, true));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read attributes", "[cooler][short]") {
  auto path = datadir / "cooler_test_file.cool";
  auto f = File::open_read_only(path.string());

  SECTION("bin size") { CHECK(f.bin_size() == 100'000); }

  SECTION("common attributes") {
    const auto& attrs = f.attributes();
    CHECK(attrs.bin_size == 100'000);
    CHECK(attrs.bin_type == "fixed");
    CHECK(attrs.creation_date == "2020-07-08T13:41:20.376258");
    CHECK(attrs.format == COOL_MAGIC);
    CHECK(attrs.format_url == "https://github.com/mirnylab/cooler");
    CHECK(attrs.format_version == 3);
    CHECK(attrs.generated_by == "cooler-0.8.8-dev");
    CHECK(attrs.assembly == "unknown");
    CHECK(attrs.metadata == "{}");
    CHECK(attrs.nbins == 26398);
    CHECK(attrs.nchroms == 20);
    CHECK(attrs.nnz == 107041);
    CHECK(attrs.storage_mode == "symmetric-upper");
    CHECK(attrs.sum.has_value());
    if (attrs.sum.has_value()) {
      std::visit([](auto& sum) { CHECK(sum == 395465); }, *attrs.sum);
    }
    CHECK(!attrs.cis.has_value());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read/write chromosomes", "[cooler][short]") {
  const auto path = (testdir() / "test_write_chroms.cool").string();

  constexpr std::uint32_t bin_size = 5000;
  const Reference chroms{Chromosome{0, "chr1", 50001}, Chromosome{1, "chr2", 25017},
                         Chromosome{2, "chr3", 10000}};

  {
    auto f = File::create_new_cooler(path, chroms, bin_size, true);
    CHECK(chroms == f.chromosomes());
  }

  const auto f = File::open_read_only(path, DEFAULT_HDF5_CACHE_SIZE, false);
  CHECK(chroms == f.chromosomes());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read/write bin table", "[cooler][short]") {
  const auto path = (testdir() / "test_write_bin_table.cool").string();

  const Reference chroms{Chromosome{0, "chr1", 50001}, Chromosome{1, "chr2", 25017},
                         Chromosome{2, "chr3", 10000}};

  constexpr std::uint32_t bin_size = 5000;
  const BinTable table(chroms, bin_size);

  { auto f = File::create_new_cooler(path, chroms, bin_size, true); }

  auto f = File::open_read_only(path);

  auto start_it = f.dataset("bins/start").begin<std::uint32_t>();
  auto end_it = f.dataset("bins/end").begin<std::uint32_t>();

  REQUIRE(start_it != f.dataset("bins/start").end<std::uint32_t>());
  REQUIRE(end_it != f.dataset("bins/end").end<std::uint32_t>());

  for (const auto bin : table) {
    CHECK(*start_it++ == bin.start());
    CHECK(*end_it++ == bin.end());
  }

  CHECK(start_it == f.dataset("bins/start").end<std::uint32_t>());
  CHECK(end_it == f.dataset("bins/end").end<std::uint32_t>());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read/write pixels", "[cooler][long]") {
  auto path1 = datadir / "cooler_test_file.cool";
  auto path2 = testdir() / "cooler_test_read_write_pixels.cool";

  using T = std::int32_t;
  auto f1 = File::open_read_only(path1.string());
  {
    auto f2 = File::create_new_cooler<T>(path2.string(), f1.chromosomes(), f1.bin_size(), true);

    const std::vector<Pixel<T>> expected(f1.begin<T>(), f1.end<T>());
    REQUIRE(expected.size() == 107041);

    std::random_device rd;
    std::mt19937_64 rand_eng{rd()};

    auto pixel_it = expected.begin();
    do {
      const auto diff = std::distance(pixel_it, expected.end());
      // Write pixels in chunks of random size
      const auto offset =
          (std::min)(diff, std::uniform_int_distribution<std::ptrdiff_t>{500, 5000}(rand_eng));
      // fmt::print(stderr, FMT_STRING("Processing {}-{} out of {}\n"),
      //            std::distance(expected.begin(), pixel_it),
      //            std::distance(expected.begin(), pixel_it + offset), expected.size());

      f2.append_pixels(pixel_it, pixel_it + offset, true);
      pixel_it += offset;
    } while (pixel_it != expected.end());
  }

  auto f2 = File::open_read_only(path2.string());

  SECTION("compare chromosomes") { CHECK(f1.chromosomes() == f2.chromosomes()); }

  SECTION("compare bins") { CHECK(f1.bins() == f2.bins()); }

  SECTION("compare indexes") {
    {
      const auto expected_chrom_offset =
          f1.dataset("indexes/chrom_offset").read_all<std::vector<std::uint64_t>>();
      const auto chrom_offset =
          f2.dataset("indexes/chrom_offset").read_all<std::vector<std::uint64_t>>();
      REQUIRE(chrom_offset.size() == expected_chrom_offset.size());
      for (std::size_t i = 0; i < chrom_offset.size(); ++i) {
        CHECK(chrom_offset[i] == expected_chrom_offset[i]);
      }
    }
    const auto expected_bin1_offset =
        f1.dataset("indexes/bin1_offset").read_all<std::vector<std::uint64_t>>();
    const auto bin1_offset =
        f2.dataset("indexes/bin1_offset").read_all<std::vector<std::uint64_t>>();
    REQUIRE(bin1_offset.size() == expected_bin1_offset.size());
    for (std::size_t i = 0; i < bin1_offset.size(); ++i) {
      CHECK(bin1_offset[i] == expected_bin1_offset[i]);
    }
  }

  SECTION("compare pixels") {
    const std::vector<Pixel<T>> expected_pixels(f1.begin<T>(), f1.end<T>());
    const std::vector<Pixel<T>> pixels(f2.begin<T>(), f2.end<T>());

    REQUIRE(expected_pixels.size() == pixels.size());
    for (std::size_t i = 0; i < pixels.size(); ++i) {
      CHECK(pixels[i] == expected_pixels[i]);
    }
  }

  SECTION("compare attributes") {
    CHECK(f1.attributes().bin_size == f2.attributes().bin_size);
    CHECK(f1.attributes().bin_type == f2.attributes().bin_type);
    CHECK(f1.attributes().format == f2.attributes().format);
    CHECK(f1.attributes().storage_mode == f2.attributes().storage_mode);
    CHECK(f1.attributes().creation_date != f2.attributes().creation_date);
    CHECK(f1.attributes().generated_by != f2.attributes().generated_by);
    CHECK(f1.attributes().assembly == f2.attributes().assembly);
    CHECK(f2.attributes().metadata == "{}");
    // Test file is still using https://github.com/mirnylab/cooler as format_url
    // CHECK(f1.attributes().format_url == f2.attributes().format_url);
    CHECK(f1.attributes().nbins == f2.attributes().nbins);
    CHECK(f1.attributes().nnz == f2.attributes().nnz);
    CHECK(f1.attributes().sum == f2.attributes().sum);
    CHECK(f2.attributes().cis == StandardAttributes::SumVar(std::int64_t(329276)));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: write weights", "[cooler][short]") {
  auto path1 = datadir / "cooler_test_file.cool";
  auto path2 = testdir() / "cooler_test_write_weights1.cool";
  auto path3 = testdir() / "cooler_test_write_weights2.cool";

  std::filesystem::remove(path2);
  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path2);
  REQUIRE_FALSE(File::open_read_only(path2.string()).has_weights("weight"));

  const auto num_bins = File::open_read_only(path1.string()).bins().size();

  SECTION("correct shape") {
    const std::vector<double> weights(num_bins, 1.23);
    File::write_weights(path2.string(), "weight", weights.begin(), weights.end());

    const auto w = *File::open_read_only(path2.string()).read_weights("weight");
    CHECK(w().size() == weights.size());
  }

  SECTION("incorrect shape") {
    std::vector<double> weights{};
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));

    weights.resize(num_bins - 1);
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));

    weights.resize(num_bins + 1);
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));
  }

  SECTION("invalid name") {
    std::vector<double> weights{};
    CHECK_THROWS(File::write_weights(path2.string(), "", weights.begin(), weights.end()));
  }

  SECTION("overwriting") {
    const std::vector<double> weights(num_bins, 1.23);
    File::write_weights(path2.string(), "weight", weights.begin(), weights.end());

    File::write_weights(path2.string(), "weight", weights.begin(), weights.end(), true);

    CHECK_THROWS(
        File::write_weights(path2.string(), "weight", weights.begin(), weights.end(), false));
  }

  SECTION("write on file creation") {
    const auto fin = File::open_read_only(path1.string());
    auto fout = File::create_new_cooler(path3.string(), fin.chromosomes(), fin.bin_size());

    const std::vector<double> weights(num_bins, 1.23);
    fout.write_weights("weight", weights.begin(), weights.end());
    fout.write_weights("weight2", weights.begin(), weights.end());
  }

  SECTION("attempt write on read-only file") {
    constexpr std::array<double, 1> w{};
    CHECK_THROWS(File::open_read_only(path2.string()).write_weights("weights", w.begin(), w.end()));
  }
}

}  // namespace hictk::cooler::test::cooler_file
