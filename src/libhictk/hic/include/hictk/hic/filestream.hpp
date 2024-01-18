// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <cstddef>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace hictk::hic::internal::filestream {

class FileStream {
  std::string _path{};
  mutable std::ifstream _ifs{};
  mutable std::ofstream _ofs{};
  std::size_t _file_size{};

 public:
  FileStream() = default;
  explicit FileStream(std::string path, std::ios::openmode mode = std::ios::in);
  static FileStream create(std::string path);

  [[nodiscard]] const std::string &path() const noexcept;
  [[nodiscard]] const std::string &url() const noexcept;
  [[nodiscard]] std::size_t size() const;

  void seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  [[nodiscard]] std::size_t tellg() const noexcept;

  void seekp(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  [[nodiscard]] std::size_t tellp() const noexcept;

  [[nodiscard]] bool eof() const noexcept;

  void flush();

  void read(std::string &buffer, std::size_t count);
  void read(char *buffer, std::size_t count);
  void read_append(std::string &buffer, std::size_t count);

  bool getline(std::string &buffer, char delim = '\n');
  [[nodiscard]] std::string getline(char delim = '\n');

  void write(std::string_view buffer);
  void write(const char *buffer, std::size_t count);

  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  [[nodiscard]] T read();
  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void read(T &buffer);

  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void write(T buffer);

  template <typename Tin,
            typename Tout = std::make_signed_t<Tin>,  // NOLINTNEXTLINE(modernize-type-traits)
            typename std::enable_if<std::is_integral<Tin>::value>::type * = nullptr>
  [[nodiscard]] Tout read_as_signed();
  template <typename Tin,
            typename Tout = std::make_unsigned_t<Tin>,  // NOLINTNEXTLINE(modernize-type-traits)
            typename std::enable_if<std::is_integral<Tin>::value>::type * = nullptr>
  [[nodiscard]] Tout read_as_unsigned();
  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename Tin, typename std::enable_if<std::is_arithmetic<Tin>::value>::type * = nullptr>
  [[nodiscard]] double read_as_double();

  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void read(std::vector<T> &buffer);

  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void write(const std::vector<T> &buffer);

  // NOLINTNEXTLINE(modernize-type-traits)
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  [[nodiscard]] std::vector<T> read(std::size_t size);

 private:
  [[nodiscard]] std::streampos new_posg(std::streamoff offset, std::ios::seekdir way);
  [[nodiscard]] std::streampos new_posp(std::streamoff offset, std::ios::seekdir way);
  void update_file_size();
  [[nodiscard]] static std::ifstream open_file_read(const std::string &path,
                                                    std::ifstream::openmode mode);
  [[nodiscard]] static std::ofstream open_file_write(const std::string &path,
                                                     std::ofstream::openmode mode);
};
}  // namespace hictk::hic::internal::filestream

#include "./impl/filestream_impl.hpp"  // NOLINT
