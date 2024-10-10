// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/pixel.hpp"
#include "hictk/transformers/pixel_merger.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

template <typename It>
using PixelMerger = ::hictk::transformers::PixelMerger<It>;

struct Coords {
  std::uint64_t bin1{};  // NOLINT
  std::uint64_t bin2{};  // NOLINT

  bool operator==(const Coords& other) const noexcept {
    return bin1 == other.bin1 && bin2 == other.bin2;
  }
  bool operator<(const Coords& other) const noexcept {
    if (bin1 == other.bin1) {
      return bin2 < other.bin2;
    }
    return bin1 < other.bin1;
  }
};

template <typename PixelIt>
static phmap::btree_map<Coords, std::int32_t> merge_pixels_hashmap(
    const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails) {
  assert(heads.size() == tails.size());

  phmap::btree_map<Coords, std::int32_t> map{};
  for (std::size_t i = 0; i < heads.size(); ++i) {
    std::for_each(heads[i], tails[i], [&](const ThinPixel<std::int32_t>& p) {
      const Coords c{p.bin1_id, p.bin2_id};

      if (map.contains(c)) {
        map[c] += p.count;
      } else {
        map[c] = p.count;
      }
    });
  }
  return map;
}

TEST_CASE("Transformers (cooler): pixel merger", "[transformers][short]") {
  const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";

  const cooler::File clr(path.string());
  const auto sel1 = clr.fetch("chr1:0-100,000,000");
  const auto sel2 = clr.fetch("chr1:50,000,000-150,000,000");
  const auto sel3 = clr.fetch("chr2:50,000,000-150,000,000");

  using It = decltype(sel1.template begin<std::int32_t>());

  SECTION("range with data") {
    const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                                sel2.template begin<std::int32_t>(),
                                sel3.template begin<std::int32_t>()};
    const std::vector<It> tails{sel1.template end<std::int32_t>(),
                                sel2.template end<std::int32_t>(),
                                sel3.template end<std::int32_t>()};

    const PixelMerger<It> merger(heads, tails);
    const auto pixels = PixelMerger<It>(heads, tails).read_all();
    const auto expected_pixels = merge_pixels_hashmap(heads, tails);

    REQUIRE(pixels.size() == expected_pixels.size());

    for (const auto& p : pixels) {
      REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
      CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
    }
  }

  SECTION("one iterator") {
    const std::vector<It> heads{sel1.template begin<std::int32_t>()};
    const std::vector<It> tails{sel1.template end<std::int32_t>()};
    const PixelMerger<It> merger(heads, tails);
    const auto pixels = PixelMerger<It>(heads, tails).read_all();
    const auto expected_pixels = merge_pixels_hashmap(heads, tails);

    REQUIRE(pixels.size() == expected_pixels.size());
    for (const auto& p : pixels) {
      REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
      CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
    }
  }

  SECTION("empty range") {
    const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                                sel2.template end<std::int32_t>(),
                                sel3.template begin<std::int32_t>()};
    const std::vector<It> tails{sel1.template end<std::int32_t>(),
                                sel2.template end<std::int32_t>(),
                                sel3.template end<std::int32_t>()};
    const PixelMerger<It> merger(heads, tails);
    const auto pixels = PixelMerger<It>(heads, tails).read_all();
    const auto expected_pixels = merge_pixels_hashmap(heads, tails);

    REQUIRE(pixels.size() == expected_pixels.size());
    for (const auto& p : pixels) {
      REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
      CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
    }
  }

  SECTION("no iterators") {
    const std::vector<It> heads{};
    const std::vector<It> tails{};
    const PixelMerger<It> merger(heads, tails);

    CHECK(merger.begin() == merger.end());
  }
}

TEST_CASE("Transformers (hic): pixel merger", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";

  const hic::File hf(path.string(), 100'000);
  const auto sel1 = hf.fetch("chr2L:0-10,000,000");
  const auto sel2 = hf.fetch("chr2L:5,000,000-15,000,000");
  const auto sel3 = hf.fetch("chr2R:5,000,000-15,000,000");

  using It = decltype(sel1.template begin<std::int32_t>());
  const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                              sel2.template begin<std::int32_t>(),
                              sel3.template begin<std::int32_t>()};
  const std::vector<It> tails{sel1.template end<std::int32_t>(), sel2.template end<std::int32_t>(),
                              sel3.template end<std::int32_t>()};

  const PixelMerger<It> merger(heads, tails);
  const auto pixels = PixelMerger<It>(heads, tails).read_all();
  const auto expected_pixels = merge_pixels_hashmap(heads, tails);

  REQUIRE(pixels.size() == expected_pixels.size());

  for (const auto& p : pixels) {
    REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
    CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::transformers
