// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/filestream.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "catch2/catch_test_macros.hpp"
#include "hictk/suppress_warnings.hpp"

using namespace hictk::hic::internal::filestream;

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

const auto path_plaintext = (hictk::test::datadir / "data.txt").string();  // NOLINT(cert-err58-cpp)
const auto path_binary = (hictk::test::datadir / "data.zip").string();     // NOLINT(cert-err58-cpp)
const auto& path = path_plaintext;

static std::string read_file(const std::string& path_) {
  std::ifstream ifs(path_, std::ios::ate);
  REQUIRE(ifs);
  const auto size = ifs.tellg();
  ifs.seekg(std::ios::beg);
  std::string buff(static_cast<std::size_t>(size), '\0');
  ifs.read(&buff.front(), size);

  return buff;
}

static std::vector<std::string> read_file_by_line(const std::string& path_, char delim = '\n') {
  std::ifstream ifs(path_);
  REQUIRE(ifs);

  std::vector<std::string> lines;
  std::string line;

  while (std::getline(ifs, line, delim)) {
    lines.push_back(line);
  }

  REQUIRE(!lines.empty());
  return lines;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream ctor", "[hic][short]") {
  SECTION("default") {
    const FileStream s{};
    CHECK(s.url().empty());
    CHECK(s.size() == 0);
  }

  SECTION("valid path") {
    const FileStream s(path_plaintext);
    CHECK(s.url() == path_plaintext);
    CHECK(s.size() == 502941);
    CHECK(!s.eof());
  }

  SECTION("invalid path") { CHECK_THROWS(FileStream("not-a-path")); }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream seek", "[hic][short]") {
  FileStream s(path_plaintext);
  {
    std::string buff;
    s.read(buff, 1);
    s.seekg(0);
  }
  SECTION("seek within chunk") {
    s.seekg(5);
    CHECK(s.tellg() == 5);

    s.seekg(10);
    CHECK(s.tellg() == 10);
  }

  SECTION("negative seek from beg") { CHECK_THROWS(s.seekg(-10)); }

  SECTION("seek from current") {
    s.seekg(10);
    CHECK(s.tellg() == 10);

    s.seekg(10, std::ios::cur);
    CHECK(s.tellg() == 20);

    s.seekg(-10, std::ios::cur);
    CHECK(s.tellg() == 10);
  }

  SECTION("seek at end") {
    s.seekg(0, std::ios::end);
    CHECK(!s.eof());
  }

  SECTION("seek past end") {
    s.seekg(0, std::ios::end);
    CHECK_THROWS(s.seekg(1, std::ios::cur));

    s.seekg(0);
    CHECK_THROWS(s.seekg(1, std::ios::end));
    CHECK(s.tellg() == 0);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream read", "[hic][short]") {
  FileStream s(path_plaintext);

  std::string buffer{"garbage"};
  const auto expected = read_file(path);
  REQUIRE(s.size() == expected.size());

  SECTION("small read") {
    s.read(buffer, 10);
    CHECK(buffer == expected.substr(0, 10));
  }

  SECTION("large read") {
    s.read(buffer, s.size());
    CHECK(buffer == expected);
  }

  SECTION("no-op read") {
    s.read(buffer, 0);
    CHECK(buffer.empty());
  }

  SECTION("seek and read") {
    const auto offset = s.size() - 10;
    s.seekg(std::int64_t(offset));
    s.read(buffer, 10);
    CHECK(buffer == expected.substr(offset));
  }

  SECTION("seek and read out-of-bound") {
    const auto offset = s.size() - 10;
    s.seekg(std::int64_t(offset));
    CHECK_THROWS(s.read(buffer, 11));
  }

  SECTION("read within chunk") {
    s.seekg(0);
    s.read(buffer, 5);
    CHECK(buffer == expected.substr(0, 5));
    s.read(buffer, 5);
    CHECK(buffer == expected.substr(5, 5));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream append", "[hic][short]") {
  FileStream s(path_plaintext);

  std::string buffer;
  const auto expected = read_file(path);

  SECTION("append to empty buffer") {
    s.append(buffer, 10);
    CHECK(buffer == expected.substr(0, 10));
  }

  SECTION("append to dirty buffer") {
    buffer = "garbage";
    s.append(buffer, 10);
    CHECK(buffer == "garbage" + expected.substr(0, 10));
  }

  SECTION("large append") {
    s.append(buffer, s.size());
    CHECK(buffer == expected);
  }

  SECTION("no-op append") {
    s.append(buffer, 0);
    CHECK(buffer.empty());
  }

  SECTION("out-of-bound read") {
    s.seekg(-1, std::ios::end);
    CHECK_THROWS(s.append(buffer, 10));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream getline", "[hic][short]") {
  FileStream s(path_plaintext);

  std::string buffer;
  const auto expected = read_file_by_line(path);

  SECTION("get one line") {
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == expected.at(0));
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == expected.at(1));
  }

  SECTION("seek and getline") {
    s.seekg(765);
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == "ibes the overall architecture of HTTP,");

    s.seekg(0);
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == expected.at(0));
  }

  SECTION("get all lines") {
    for (std::size_t i = 0; s.getline(buffer); ++i) {
      CHECK(buffer == expected[i]);
    }
    CHECK(s.eof());

    CHECK_THROWS(s.getline(buffer));
    CHECK(buffer.empty());
  }

  SECTION("custom delimiter") {
    s.seekg(74);
    CHECK(s.getline(',').empty());

    CHECK(s.getline(buffer, ':'));
    CHECK(buffer == " Ed.\nRequest for Comments");
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: filestream read binary", "[hic][short]") {
  FileStream s(path_binary);
  s.seekg(10);

  // Expected numbers were computed with Python code like the following:
  // for t in [np.uint8, np.uint16, np.uint32, np.uint64, np.int8, np.int16, np.int32, np.int64,
  // np.float32, np.float64]:
  //     print(t, np.frombuffer(open("data.zip", "rb").read()[10:18], dtype=t)[0])

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_USELESS_CAST
  SECTION("uint8") { CHECK(s.read<std::uint8_t>() == std::uint8_t(162)); }
  SECTION("uint16") { CHECK(s.read<std::uint16_t>() == std::uint16_t(42658)); }
  SECTION("uint32") { CHECK(s.read<std::uint32_t>() == std::uint32_t(1433446050)); }
  SECTION("uint64") { CHECK(s.read<std::uint64_t>() == std::uint64_t(18260117889181853346ULL)); }

  SECTION("int8") { CHECK(s.read<std::int8_t>() == std::int8_t(-94)); }
  SECTION("int16") { CHECK(s.read<std::int16_t>() == std::int16_t(-22878)); }
  SECTION("int32") { CHECK(s.read<std::int32_t>() == std::int32_t(1433446050)); }
  SECTION("int64") { CHECK(s.read<std::int64_t>() == std::int64_t(-186626184527698270)); }

  SECTION("float") { CHECK(s.read<float>() == 16537405000000.0F); }
  SECTION("double") { CHECK(s.read<double>() == -1.2758357206942371e+296); }

  SECTION("char") { CHECK(s.read<char>() == static_cast<char>(162)); }
  SECTION("unsigned char") { CHECK(s.read<unsigned char>() == static_cast<unsigned char>(162)); }
  DISABLE_WARNING_POP

  SECTION("vector") {
    s.seekg(0);
    constexpr std::array<std::int32_t, 32> expected{
        67324752,    20,          -1499332600, -126266000,  316472680,   -71892991,  720898,
        926220316,   758592304,   2020879920,  156521844,   1067451136,  1101095797, 2020959093,
        67174411,    501,         5124,        -1141015552, -1772542862, 787614245,  1386282978,
        -1957338045, 1449544581,  1142046551,  -518143477,  -1249957234, 831590659,  -732484307,
        1294996684,  -1436898904, 1231094186,  1614771469};

    const auto buffer = s.read<std::int32_t>(expected.size());
    REQUIRE(expected.size() == buffer.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == buffer[i]);
    }
  }
}
