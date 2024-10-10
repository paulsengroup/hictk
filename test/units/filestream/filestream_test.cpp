// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/filestream.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "hictk/suppress_warnings.hpp"
#include "hictk/type_traits.hpp"
#include "tmpdir.hpp"

namespace hictk::filestream::test {

using namespace hictk::filestream;

const auto path_plaintext = (datadir / "data.txt").string();  // NOLINT(cert-err58-cpp)
const auto path_binary = (datadir / "data.zip").string();     // NOLINT(cert-err58-cpp)
static const auto& path = path_plaintext;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
[[maybe_unused]] static std::string read_file(const std::string& path_) {
  std::ifstream ifs(path_, std::ios::ate);
  REQUIRE(ifs);
  const auto size = ifs.tellg();
  ifs.seekg(std::ios::beg);
  std::string buff(static_cast<std::size_t>(size), '\0');
  ifs.read(&buff.front(), size);

  return buff;
}

[[maybe_unused]] static std::vector<std::string> read_file_by_line(const std::string& path_,
                                                                   char delim = '\n') {
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

TEST_CASE("FileStream ctor", "[filestream][short]") {
  SECTION("default") {
    const FileStream<> s{};
    CHECK(s.path().empty());
    CHECK(s.size() == 0);
  }

  SECTION("valid path (read)") {
    const FileStream<> s(path_plaintext, nullptr);
    CHECK(s.path() == path_plaintext);
    CHECK(s.size() == 502941);
    CHECK(!s.eof());
  }

  SECTION("valid path (write)") {
    const auto path1 = testdir() / "filestream_ctor_write.bin";
    const auto s = FileStream<>::create(path1.string(), nullptr);
    CHECK(s.path() == path1);
    CHECK(s.size() == 0);
    CHECK(!s.eof());
  }

  SECTION("invalid path") { CHECK_THROWS(FileStream<>("not-a-path", nullptr)); }
}

TEST_CASE("FileStream seek", "[filestream][short]") {
  SECTION("read") {
    FileStream<> s(path_plaintext, nullptr);
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
      CHECK_THROWS(s.seekg(-1, std::ios::end));
      CHECK(s.tellg() == 0);
    }
  }
  SECTION("write") {
    const auto path1 = testdir() / "filestream_seek.bin";
    std::filesystem::remove(path1);  // NOLINT
    auto s = FileStream<>::create(path1.string(), nullptr);

    SECTION("seek within chunk") {
      s.seekp(5);
      CHECK(s.tellp() == 5);

      s.seekp(10);
      CHECK(s.tellp() == 10);
    }

    SECTION("negative seek from beg") { CHECK_THROWS(s.seekp(-10)); }

    SECTION("seek from current") {
      s.seekp(10);
      CHECK(s.tellp() == 10);

      s.seekp(10, std::ios::cur);
      CHECK(s.tellp() == 20);

      s.seekp(-10, std::ios::cur);
      CHECK(s.tellp() == 10);
    }

    SECTION("seek at end") {
      s.seekg(0, std::ios::end);
      CHECK(!s.eof());
    }

    SECTION("seek past end") {
      s.seekp(0, std::ios::end);
      CHECK_NOTHROW(s.seekp(1, std::ios::cur));
    }
  }
}

TEST_CASE("FileStream read", "[filestream][short]") {
  FileStream<> s(path_plaintext, nullptr);
  std::string buffer{"garbage"};
  const auto expected = read_file(path);
  REQUIRE(s.size() == static_cast<std::streamsize>(expected.size()));

  SECTION("small read") {
    s.read(buffer, 10);
    CHECK(buffer == expected.substr(0, 10));
  }

  SECTION("large read") {
    s.read(buffer, static_cast<std::size_t>(s.size()));
    CHECK(buffer == expected);
  }

  SECTION("no-op read") {
    s.read(buffer, 0);
    CHECK(buffer.empty());
  }

  SECTION("seek and read") {
    const auto offset = static_cast<std::size_t>(s.size() - 10);
    s.seek_and_read(static_cast<std::int64_t>(offset), buffer, 10);
    CHECK(buffer == expected.substr(offset));
  }

  SECTION("seek and read out-of-bound") {
    const auto offset = static_cast<std::streampos>(s.size() - 10);
    CHECK_THROWS(s.seek_and_read(offset, buffer, 11));
  }

  SECTION("read within chunk") {
    s.seekg(0);
    s.read(buffer, 5);
    CHECK(buffer == expected.substr(0, 5));
    s.read(buffer, 5);
    CHECK(buffer == expected.substr(5, 5));
  }

  SECTION("multi-threaded") {
    s = FileStream<>(path_plaintext, std::make_shared<std::mutex>());

    const std::streampos offset1 = 0;
    const std::streampos offset2 = 5;

    s.seekg(offset1);
    const auto expected1 = s.read(10);
    s.seekg(offset2);
    const auto expected2 = s.read(10);

    REQUIRE(expected1.size() == 10);
    REQUIRE(expected2.size() == 10);

    std::atomic<std::size_t> threads_started{};
    std::mutex catch2_mtx{};

    auto worker = [&](std::size_t id, std::streampos offset, std::string_view expected_) {
      std::string buffer_;
      std::size_t tests = 0;
      std::size_t failures = 0;

      threads_started++;
      while (threads_started != 2);  // NOLINT

      try {
        const auto timepoint = std::chrono::steady_clock::now() + std::chrono::seconds(5);
        for (; std::chrono::steady_clock::now() < timepoint; ++tests) {
          buffer_.clear();
          const auto new_offset = s.seek_and_read(offset, buffer_, expected_.size()).second;
          // NOLINTNEXTLINE(readability-implicit-bool-conversion)
          failures += new_offset != offset + std::streamoff(expected_.size());

          [[maybe_unused]] const auto lck = std::scoped_lock(catch2_mtx);
          CHECK(new_offset == offset + std::streamoff(expected_.size()));
        }

        return std::make_pair(tests, failures);
      } catch (const std::exception& e) {
        throw std::runtime_error("Exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + "): " + e.what());
      } catch (...) {
        throw std::runtime_error("Unknown exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + ")");
      }
    };

    auto w1 = std::async(std::launch::async, worker, 1, offset1, expected1);
    auto w2 = std::async(std::launch::async, worker, 2, offset2, expected2);

    const auto [tests1, fails1] = w1.get();
    const auto [tests2, fails2] = w2.get();
    printf("performed %zu reads (%zu failures)\n", tests1 + tests2, fails1 + fails2);  // NOLINT
  }
}

TEST_CASE("FileStream read_append", "[filestream][short]") {
  FileStream<> s(path_plaintext, nullptr);

  std::string buffer;
  const auto expected = read_file(path);

  SECTION("append to empty buffer") {
    s.read_append(buffer, 10);
    CHECK(buffer == expected.substr(0, 10));
  }

  SECTION("append to dirty buffer") {
    buffer = "garbage";
    s.read_append(buffer, 10);
    CHECK(buffer == "garbage" + expected.substr(0, 10));
  }

  SECTION("large append") {
    s.read_append(buffer, static_cast<std::size_t>(s.size()));
    CHECK(buffer == expected);
  }

  SECTION("no-op append") {
    s.read_append(buffer, 0);
    CHECK(buffer.empty());
  }

  SECTION("out-of-bound read") {
    s.seekg(1, std::ios::end);
    CHECK_THROWS(s.read_append(buffer, 10));
  }
}

TEST_CASE("FileStream getline", "[filestream][short]") {
  FileStream<> s(path_plaintext, nullptr);

  std::string buffer;
  const auto expected = read_file_by_line(path);

  SECTION("get one line") {
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == expected.at(0));
    CHECK(s.getline(buffer) == true);
    CHECK(buffer == expected.at(1));
  }

  SECTION("seek and getline") {
    auto status = s.seek_and_getline(765, buffer);
    CHECK(std::get<0>(status) == true);
    CHECK(buffer == "ibes the overall architecture of HTTP,");

    status = s.seek_and_getline(0, buffer);
    CHECK(std::get<0>(status) == true);
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

  SECTION("multi-threaded") {
    s = FileStream<>(path_plaintext, std::make_shared<std::mutex>());

    const std::streampos offset1 = 25;
    const std::streampos offset2 = 30;

    s.seekg(offset1);
    const auto expected1 = s.getline();
    s.seekg(offset2);
    const auto expected2 = s.getline();

    REQUIRE(!expected1.empty());
    REQUIRE(!expected2.empty());

    std::atomic<std::size_t> threads_started{};
    std::mutex catch2_mtx{};

    auto worker = [&](std::size_t id, std::streampos offset, std::string_view expected_) {
      std::string buffer_;
      std::size_t tests = 0;
      std::size_t failures = 0;

      const auto offset_after_read_expected = offset + std::streamoff(expected_.size() + 1);

      threads_started++;
      while (threads_started != 2);  // NOLINT

      try {
        const auto timepoint = std::chrono::steady_clock::now() + std::chrono::seconds(5);
        for (; std::chrono::steady_clock::now() < timepoint; ++tests) {
          buffer_.clear();
          const auto status = s.seek_and_getline(offset, buffer_);
          const auto delimiter_found = std::get<0>(status);
          const auto offset_after_read = std::get<2>(status);

          failures += static_cast<std::size_t>(offset_after_read != offset_after_read_expected ||
                                               expected_ != buffer_);

          [[maybe_unused]] const auto lck = std::scoped_lock(catch2_mtx);
          CHECK(delimiter_found);
          CHECK(offset_after_read == offset_after_read_expected);
          CHECK(expected_ == buffer_);
        }

        return std::make_pair(tests, failures);
      } catch (const std::exception& e) {
        throw std::runtime_error("Exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + "): " + e.what());
      } catch (...) {
        throw std::runtime_error("Unknown exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + ")");
      }
    };

    auto w1 = std::async(std::launch::async, worker, 1, offset1, expected1);
    auto w2 = std::async(std::launch::async, worker, 2, offset2, expected2);

    const auto [tests1, fails1] = w1.get();
    const auto [tests2, fails2] = w2.get();
    printf("performed %zu reads (%zu failures)\n", tests1 + tests2, fails1 + fails2);  // NOLINT
  }
}

TEST_CASE("FileStream read binary", "[filestream][short]") {
  FileStream<> s(path_binary, nullptr);
  const std::streampos offset = 10;
  s.seekg(offset);

  // Expected numbers were computed with Python code like the following:
  // for t in [np.uint8, np.uint16, np.uint32, np.uint64, np.int8, np.int16, np.int32, np.int64,
  // np.float32, np.float64]:
  //     print(t, np.frombuffer(open("data.zip", "rb").read()[10:18], dtype=t)[0])

  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_USELESS_CAST
  SECTION("uint8") { CHECK(s.read<std::uint8_t>() == std::uint8_t{162}); }
  SECTION("uint16") { CHECK(s.read<std::uint16_t>() == std::uint16_t{42658}); }
  SECTION("uint32") {
    CHECK(s.read<std::uint32_t>() == std::uint32_t{1433446050});
    CHECK(std::is_same_v<decltype(s.read_as_signed<std::uint32_t>()), std::int32_t>);
    s.seekg(offset);
    CHECK(s.read_as_signed<std::uint32_t>() == std::int32_t{1433446050});
  }
  SECTION("uint64") {
    CHECK(s.read<std::uint64_t>() == std::uint64_t{18260117889181853346ULL});
    CHECK(std::is_same_v<decltype(s.read_as_signed<std::uint64_t>()), std::int64_t>);
    s.seekg(offset);
    CHECK(s.read_as_signed<std::uint64_t>() == std::int64_t(-186626184527698270LL));
  }

  SECTION("int8") { CHECK(s.read<std::int8_t>() == std::int8_t(-94)); }
  SECTION("int16") { CHECK(s.read<std::int16_t>() == std::int16_t(-22878)); }
  SECTION("int32") {
    CHECK(s.read<std::int32_t>() == std::int32_t{1433446050});
    CHECK(std::is_same_v<decltype(s.read_as_unsigned<std::int32_t>()), std::uint32_t>);
    s.seekg(offset);
    CHECK(s.read_as_unsigned<std::int32_t>() == std::uint32_t{1433446050});
  }
  SECTION("int64") {
    CHECK(s.read<std::int64_t>() == std::int64_t(-186626184527698270));
    CHECK(std::is_same_v<decltype(s.read_as_unsigned<std::int64_t>()), std::uint64_t>);
    s.seekg(offset);
    CHECK(s.read_as_unsigned<std::int64_t>() == std::uint64_t{18260117889181853346ULL});
  }

  SECTION("float") { CHECK(s.read<float>() == 16537405000000.0F); }
  SECTION("double") {
    CHECK(s.read<double>() == -1.2758357206942371e+296);
    s.seekg(offset);
    CHECK(s.read_as_double<float>() == 16537404571648.0);
  }

  SECTION("char") { CHECK(s.read<char>() == static_cast<char>(162)); }
  SECTION("unsigned char") { CHECK(s.read<unsigned char>() == static_cast<unsigned char>(162)); }
  HICTK_DISABLE_WARNING_POP

  SECTION("vector") {
    constexpr std::array<std::int32_t, 32> expected{
        67324752,    20,          -1499332600, -126266000,  316472680,   -71892991,  720898,
        926220316,   758592304,   2020879920,  156521844,   1067451136,  1101095797, 2020959093,
        67174411,    501,         5124,        -1141015552, -1772542862, 787614245,  1386282978,
        -1957338045, 1449544581,  1142046551,  -518143477,  -1249957234, 831590659,  -732484307,
        1294996684,  -1436898904, 1231094186,  1614771469};

    s.seekg(0);
    auto buffer = s.read_vector<std::int32_t>(expected.size());
    REQUIRE(expected.size() == buffer.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == buffer[i]);
    }

    std::fill(buffer.begin(), buffer.end(), 0);
    s.seek_and_read(0, buffer);
    REQUIRE(expected.size() == buffer.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == buffer[i]);
    }
  }
}

TEST_CASE("FileStream write", "[filestream][short]") {
  const auto tmpfile = testdir() / "filestream_write.bin";
  std::filesystem::remove(tmpfile);  // NOLINT
  auto s = FileStream<>::create(tmpfile.string(), nullptr);

  SECTION("small write") {
    constexpr std::string_view buffer{"test"};
    s.write(buffer);
    CHECK(s.size() == buffer.size());
  }

  SECTION("large write") {
    const auto buffer = read_file(path);
    s.write(buffer);
    CHECK(s.size() == static_cast<std::streamsize>(buffer.size()));
  }

  SECTION("no-op read") {
    constexpr std::string_view buffer{};
    s.write(buffer);
    CHECK(s.size() == 0);
  }

  SECTION("seek and write") {
    const std::size_t offset = 10;
    s.seekp(std::int64_t{offset});
    constexpr std::string_view buffer{"test"};
    s.write(buffer);
    CHECK(s.size() == buffer.size() + offset);
  }
  SECTION("multi-threaded") {
    s.close();
    std::filesystem::remove(tmpfile);  // NOLINT
    s = FileStream<>::create(tmpfile.string(), std::make_shared<std::mutex>());

    std::string_view msg1{"0123456789"};
    std::string_view msg2{"abcdefghijklmnopqrstwxyz"};

    const std::streampos offset1 = 0;
    const std::streampos offset2 = offset1 + static_cast<std::streamoff>(msg1.size());

    const std::streamsize expected_file_size = offset2 + static_cast<std::streamoff>(msg2.size());

    REQUIRE(s.size() == 0);

    std::atomic<std::size_t> threads_started{};
    std::mutex catch2_mtx{};

    auto worker = [&](std::size_t id, std::streampos offset, std::string_view message) {
      std::size_t tests = 0;
      std::size_t failures = 0;

      threads_started++;
      while (threads_started != 2);  // NOLINT

      try {
        const auto timepoint = std::chrono::steady_clock::now() + std::chrono::seconds(5);
        for (; std::chrono::steady_clock::now() < timepoint; ++tests) {
          const auto new_offset = s.seek_and_write(offset, message).second;
          // NOLINTNEXTLINE(readability-implicit-bool-conversion)
          failures += offset + static_cast<std::streamoff>(message.size()) != new_offset;

          [[maybe_unused]] const auto lck = std::scoped_lock(catch2_mtx);
          CHECK(offset + static_cast<std::streamoff>(message.size()) == new_offset);
        }

        return std::make_pair(tests, failures);
      } catch (const std::exception& e) {
        throw std::runtime_error("Exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + "): " + e.what());
      } catch (...) {
        throw std::runtime_error("Unknown exception caught in worker #" + std::to_string(id) +
                                 " (iteration " + std::to_string(tests) + ")");
      }
    };

    auto w1 = std::async(std::launch::async, worker, 1, offset1, msg1);
    auto w2 = std::async(std::launch::async, worker, 2, offset2, msg2);

    const auto [tests1, fails1] = w1.get();
    const auto [tests2, fails2] = w2.get();
    printf("performed %zu writes (%zu failures)\n", tests1 + tests2, fails1 + fails2);  // NOLINT

    REQUIRE(s.size() == expected_file_size);

    std::string buff{};
    s.seek_and_read(offset1, buff, msg1.size());
    CHECK(buff == msg1);
    s.seek_and_read(offset2, buff, msg2.size());
    CHECK(buff == msg2);
  }
}

template <typename T>
static void write_and_compare(FileStream<>& s, const T& data) {
  s.write(data);
  s.flush();
  REQUIRE(s.size() == sizeof(T));
  CHECK(s.read<T>() == data);
}

TEST_CASE("FileStream write binary", "[filestream][short]") {
  const auto tmpfile = testdir() / "filestream_write_binary.bin";
  std::filesystem::remove(tmpfile);  // NOLINT
  auto s = FileStream<>::create(tmpfile.string(), nullptr);

  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_USELESS_CAST
  SECTION("uint8") { write_and_compare(s, std::uint8_t{162}); }
  SECTION("uint16") { write_and_compare(s, std::uint16_t{42658}); }
  SECTION("uint32") { write_and_compare(s, std::uint32_t{1433446050}); }
  SECTION("uint64") { write_and_compare(s, std::uint64_t{18260117889181853346ULL}); }

  SECTION("int8") { write_and_compare(s, std::int8_t(-94)); }
  SECTION("int16") { write_and_compare(s, std::int16_t(-22878)); }
  SECTION("int32") { write_and_compare(s, std::int32_t{1433446050}); }
  SECTION("int64") { write_and_compare(s, std::int64_t(-186626184527698270)); }

  SECTION("float") { write_and_compare(s, 16537405000000.0F); }
  SECTION("double") { write_and_compare(s, -1.2758357206942371e+296); }

  SECTION("bool") { write_and_compare(s, false); }
  SECTION("char") { write_and_compare(s, static_cast<char>(162)); }
  SECTION("unsigned char") { write_and_compare(s, static_cast<unsigned char>(162)); }
  HICTK_DISABLE_WARNING_POP

  SECTION("vector") {
    std::vector<std::int32_t> data{
        67324752,    20,          -1499332600, -126266000,  316472680,   -71892991,  720898,
        926220316,   758592304,   2020879920,  156521844,   1067451136,  1101095797, 2020959093,
        67174411,    501,         5124,        -1141015552, -1772542862, 787614245,  1386282978,
        -1957338045, 1449544581,  1142046551,  -518143477,  -1249957234, 831590659,  -732484307,
        1294996684,  -1436898904, 1231094186,  1614771469};

    s.write(data);
    s.flush();
    REQUIRE(s.size() == static_cast<std::streamsize>(sizeof(std::int32_t) * data.size()));
    const auto buffer = s.read_vector<std::int32_t>(data.size());
    REQUIRE(data.size() == buffer.size());
    for (std::size_t i = 0; i < data.size(); ++i) {
      CHECK(data[i] == buffer[i]);
    }
  }
}

TEST_CASE("FileStream resize", "[filestream][short]") {
  const auto tmpfile = testdir() / "filestream_write.bin";
  std::filesystem::remove(tmpfile);  // NOLINT
  auto s = FileStream<>::create(tmpfile.string(), nullptr);

  const std::string_view msg{"this is a relatively long string"};

  s.write(msg);
  CHECK(s.size() == conditional_static_cast<std::streamsize>(msg.size()));
  CHECK(s.tellg() == 0);
  CHECK(s.tellp() == s.size());

  s.resize(5);
  CHECK(s.size() == std::streamsize{5});
  CHECK(s.tellg() == 0);
  CHECK(s.tellp() == 5);

  s.resize(100);
  CHECK(s.size() == std::streamsize{100});
  CHECK(s.tellg() == 0);
  CHECK(s.tellp() == 5);
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::filestream::test
