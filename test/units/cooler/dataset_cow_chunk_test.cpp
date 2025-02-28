// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include "hictk/cooler/dataset.hpp"

namespace hictk::cooler::test::dataset {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
template <typename T>
static void validate_chunk(const internal::COWChunk<T>& chunk, const std::vector<T>& data) {
  REQUIRE(chunk.size() == data.size());

  if (chunk.empty()) {
    return;
  }

  const auto& buff = chunk();

  for (std::size_t i = 0; i < chunk.size(); ++i) {
    if constexpr (std::is_floating_point_v<T>) {
      CHECK_THAT(buff[i], Catch::Matchers::WithinRel(data[i]));
    } else {
      CHECK(buff[i] == data[i]);
    }
  }
}

[[nodiscard]] static bool shrink_to_fit_is_honored() {
  std::vector<std::uint64_t> v(10);
  v.resize(5);
  return v.capacity() == v.size();
}

TEST_CASE("Cooler: dataset COWChunk", "[dataset][short]") {
  using COWChunk = internal::COWChunk<std::uint64_t>;
  SECTION("ctors") {
    SECTION("default") {
      const COWChunk chunk{};
      CHECK(chunk.empty());
      CHECK(chunk.capacity() == 0);

      CHECK(chunk.id() == 0);
      CHECK(chunk.start() == 0);
      CHECK(chunk.end() == 0);
    }

    const std::vector<std::uint64_t> vec{1, 2, 3};
    SECTION("vector") {
      const COWChunk chunk{10, vec, 5};

      CHECK(chunk.size() == 3);
      CHECK(chunk.capacity() == 5);

      CHECK(chunk.id() == 2);
      CHECK(chunk.start() == 10);
      CHECK(chunk.end() == 13);

      CHECK(COWChunk{10, vec}.capacity() == vec.size());
    }
  }

  SECTION("copy on write") {
    const std::vector<std::uint64_t> vec{1, 2, 3};

    SECTION("default constructed") {
      const COWChunk chunk{};
      CHECK(chunk.use_count() == 0);
    }

    SECTION("copy assignment") {
      const COWChunk chunk1{0, vec};
      CHECK(chunk1.use_count() == 1);
      {
        const auto chunk2 = chunk1;  // NOLINT
        CHECK(chunk1.use_count() == 2);
        CHECK(&chunk1() == &chunk2());
      }
      CHECK(chunk1.use_count() == 1);
    }

    SECTION("move assignment") {
      COWChunk chunk1{0, vec};
      CHECK(chunk1.use_count() == 1);

      const auto chunk2 = std::move(chunk1);
      CHECK(chunk2.use_count() == 1);
      CHECK(chunk1.use_count() == 0);  // NOLINT(*-use-after-move, *-invalid-access-moved)
    }
  }

  SECTION("accessors") {
    const std::vector<std::uint64_t> vec{1, 2, 3};
    const COWChunk chunk{10, vec, 5};

    for (std::size_t i = 0; i < 20; ++i) {
      if (i < chunk.start() || i >= chunk.end()) {
        CHECK(!chunk(i).has_value());

      } else {
        CHECK(chunk(i) == vec[i - chunk.start()]);
        CHECK(chunk[i] == vec[i - chunk.start()]);
      }
    }
  }

  SECTION("update") {
    const std::vector<std::uint64_t> vec{1, 2, 3};
    COWChunk chunk{0, vec};

    SECTION("start only") {
      REQUIRE(chunk.start() == 0);
      REQUIRE(chunk.end() == 3);
      REQUIRE(chunk.id() == 0);
      chunk.update(10);
      CHECK(chunk.start() == 10);
      CHECK(chunk.end() == 13);
      CHECK(chunk.id() == 3);
    }

    SECTION("shrink") {
      const std::vector<std::uint64_t> new_vec{10, 20};
      const auto* old_buff = &chunk();

      REQUIRE(chunk.size() == 3);
      REQUIRE(chunk.capacity() == 3);
      chunk.update(0, new_vec);

      CHECK(chunk.size() == 2);
      CHECK(chunk.capacity() == 3);
      CHECK(&chunk() == old_buff);
      validate_chunk(chunk, new_vec);
    }

    SECTION("no size change") {
      const std::vector<std::uint64_t> new_vec{10, 20, 30};
      const auto* old_buff = &chunk();

      REQUIRE(chunk.size() == 3);
      REQUIRE(chunk.capacity() == 3);
      chunk.update(0, new_vec);

      CHECK(chunk.size() == 3);
      CHECK(chunk.capacity() == 3);
      CHECK(&chunk() == old_buff);
      validate_chunk(chunk, new_vec);
    }

    SECTION("grow") {
      const std::vector<std::uint64_t> new_vec{10, 20, 30, 40};

      REQUIRE(chunk.size() == 3);
      REQUIRE(chunk.capacity() == 3);

      CHECK_THROWS_WITH(
          chunk.update(0, new_vec),
          Catch::Matchers::ContainsSubstring("incoming data is larger than the available space"));
      CHECK_THROWS_WITH(chunk.update(0, std::make_shared<std::vector<std::uint64_t>>(new_vec)),
                        Catch::Matchers::ContainsSubstring(
                            "incoming data has a different size then the current buffer"));
      chunk.reserve(new_vec.size());
      chunk.update(0, new_vec);

      CHECK(chunk.size() == 4);
      CHECK(chunk.capacity() == 4);
      validate_chunk(chunk, new_vec);
    }

    SECTION("empty vector") {
      REQUIRE(chunk.size() == 3);
      REQUIRE(chunk.capacity() == 3);
      chunk.update(0, std::vector<std::uint64_t>{});

      CHECK(chunk.empty());
      CHECK(chunk.capacity() == 0);
      // only the const overload of operator() can be called on an empty buffer
      const auto const_chunk = chunk;
      CHECK(const_chunk().empty());
    }

    SECTION("reallocation required") {
      const auto vec_ptr = std::make_shared<std::vector<std::uint64_t>>(vec);
      const std::vector<std::uint64_t> new_vec{10, 20};
      chunk = COWChunk{0, vec_ptr};
      const auto* old_buff = &chunk();

      REQUIRE(vec_ptr.use_count() == 2);
      REQUIRE(chunk.use_count() == 2);

      chunk.update(0, new_vec);

      CHECK(chunk.size() == 2);
      CHECK(chunk.capacity() == 3);
      CHECK(chunk.use_count() == 1);
      CHECK(&chunk() != old_buff);
      validate_chunk(chunk, new_vec);
    }
  }

  SECTION("resize") {
    const std::vector<std::uint64_t> vec{1, 2, 3};
    COWChunk chunk{0, vec};
    const auto* old_buff = &chunk();

    REQUIRE(chunk.size() == 3);
    REQUIRE(chunk.capacity() == 3);

    SECTION("no-op") {
      chunk.resize(3);
      CHECK(chunk.size() == 3);
      CHECK(chunk.capacity() == 3);
      CHECK(old_buff == &chunk());
    }

    SECTION("shrink") {
      chunk.resize(1);
      CHECK(chunk.size() == 1);
      CHECK(chunk.capacity() == 3);
      CHECK(old_buff == &chunk());
    }

    SECTION("shrink to zero") {
      chunk.resize(0);
      CHECK(chunk.empty());
      CHECK(chunk.capacity() == 0);

      CHECK(chunk().empty());
    }

    SECTION("grow") {
      chunk.resize(10);
      CHECK(chunk.size() == 10);
      CHECK(chunk.capacity() == 10);
    }

    SECTION("reallocation required") {
      const auto vec_ptr = std::make_shared<std::vector<std::uint64_t>>(vec);
      chunk = COWChunk{0, vec_ptr};

      REQUIRE(vec_ptr.use_count() == 2);
      REQUIRE(chunk.use_count() == 2);

      chunk.resize(10);

      CHECK(chunk.size() == 10);
      CHECK(chunk.capacity() == 10);
      CHECK(chunk.use_count() == 1);

      chunk.resize(5);

      CHECK(chunk.size() == 5);
      CHECK(chunk.capacity() == 10);
      CHECK(chunk.use_count() == 1);
      if (chunk.empty()) {
        for (std::size_t i = 0; i < std::min(chunk.size(), vec.size()); ++i) {
          CHECK(chunk[i] == vec[i]);
        }
      }

      if (shrink_to_fit_is_honored()) {
        chunk.resize(5, true);
        CHECK(chunk.size() == 5);
        CHECK(chunk.capacity() == 5);
        CHECK(chunk.use_count() == 1);
      }
    }
  }

  SECTION("reserve") {
    COWChunk chunk{0, std::vector<std::uint64_t>{}};

    chunk.reserve(10);
    REQUIRE(chunk.capacity() == 10);
    chunk.reserve(5);
    REQUIRE(chunk.capacity() == 10);
    chunk.reserve(20);
    REQUIRE(chunk.capacity() == 20);
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::dataset
