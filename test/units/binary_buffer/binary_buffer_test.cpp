// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/binary_buffer.hpp"

#include <catch2/catch_test_macros.hpp>

namespace hictk::binary_buffer::test {

TEST_CASE("BinaryBuffer") {
  const std::vector<std::uint32_t> ibuff{1, 2, 3};
  // NOLINTNEXTLINE
  const std::string sibuff{reinterpret_cast<const char*>(ibuff.data()),
                           ibuff.size() * sizeof(std::uint32_t)};
  const std::string sbuff{"Hi\nThere\0!"};  // NOLINT

  SECTION("read") {
    BinaryBuffer buff{};

    SECTION("std::string") {
      buff.reset() = sbuff;
      std::string read_buff(2, '\0');

      buff.read(read_buff, read_buff.size());
      CHECK(read_buff == "Hi");

      buff.read(read_buff, 1);
      CHECK(read_buff == "\n");
    }

    SECTION("char*") {
      buff.reset() = sbuff;
      std::string read_buff(2, '\0');

      buff.read(read_buff.data(), read_buff.size());
      CHECK(read_buff == "Hi");

      buff.read(read_buff, 1);
      CHECK(read_buff == "\n");
    }

    SECTION("N") {
      buff.reset() = sibuff;
      CHECK(buff.read<std::uint32_t>() == 1);
      CHECK(buff.read<std::uint32_t>() == 2);
      CHECK(buff.read<std::uint32_t>() == 3);
    }

    SECTION("std::vector<N>") {
      std::vector<std::uint32_t> read_buff(3, 0);

      buff.reset() = sibuff;
      buff.read(read_buff);
      REQUIRE(read_buff.size() == 3);
      CHECK(read_buff[0] == 1);
      CHECK(read_buff[1] == 2);
      CHECK(read_buff[2] == 3);
    }

    SECTION("getline") {
      SECTION("newline") {
        buff.reset() = sbuff;
        CHECK(buff.getline() == "Hi");
      }
      SECTION("nullterm") {
        buff.reset() = sbuff;
        CHECK(buff.getline('\0') == "Hi\nThere");
      }
    }
  }

  SECTION("write") {
    BinaryBuffer buff{};

    SECTION("std::string") {
      buff.write("test");
      std::string read_buff(4, '\0');
      buff.read(read_buff, read_buff.size());
      CHECK(read_buff == "test");
    }

    SECTION("N") {
      buff.write(std::int64_t(123));
      CHECK(buff.read<std::int64_t>() == 123);
    }

    SECTION("std::vector<N>") {
      buff.write(ibuff);
      std::vector<std::uint32_t> read_buff(3, 0);
      buff.read(read_buff);

      REQUIRE(ibuff.size() == read_buff.size());
      CHECK(ibuff[0] == read_buff[0]);
      CHECK(ibuff[1] == read_buff[1]);
      CHECK(ibuff[2] == read_buff[2]);
    }
  }
}

}  // namespace hictk::binary_buffer::test
