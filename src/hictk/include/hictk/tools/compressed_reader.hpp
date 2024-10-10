// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <archive.h>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

struct archive;
struct archive_entry;

namespace hictk::tools::io {

class CompressedReader {
  using archive_ptr_t = std::unique_ptr<archive, decltype(&archive_read_free)>;

 public:
  CompressedReader() = default;
  explicit CompressedReader(
      const std::filesystem::path& path,  // NOLINTNEXTLINE(*-avoid-magic-numbers)
      std::size_t buff_capacity = 512UL << 10U);

  bool getline(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string_view getline(char sep = '\n');
  bool readall(std::string& buff, char sep = '\n');
  [[nodiscard]] std::string readall(char sep = '\n');
  [[nodiscard]] bool eof() const noexcept;
  [[nodiscard]] bool is_open() const noexcept;
  void close();
  void open(const std::filesystem::path& path);
  void reset();

  [[nodiscard]] explicit operator bool() const;
  [[nodiscard]] bool operator!() const;

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
  [[nodiscard]] std::string path_string() const noexcept;

 private:
  std::filesystem::path _path{};
  archive_ptr_t _arc{nullptr, archive_read_free};
  std::unique_ptr<archive_entry*> _arc_entry{new archive_entry* };
  std::string _buff{};
  std::string _tok_tmp_buff{};
  std::size_t _idx{0};
  bool _eof{false};

  void handle_libarchive_errors(la_ssize_t errcode) const;
  void handle_libarchive_errors() const;
  // Returns false when reaching eof
  [[nodiscard]] bool read_next_chunk();
  // Returns false when unable to find the next token occurrence
  [[nodiscard]] bool read_next_token(std::string& buff, char sep);
  [[nodiscard]] std::string_view read_next_token(char sep);
};

}  // namespace hictk::tools::io
