// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <random>

#include "hictk/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

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
          std::min(diff, std::uniform_int_distribution<std::ptrdiff_t>{500, 5000}(rand_eng));
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

}  // namespace hictk::cooler::test::cooler_file
