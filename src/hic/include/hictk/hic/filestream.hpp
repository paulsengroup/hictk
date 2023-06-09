// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictk/hic/common.hpp"

namespace hictk::internal::filestream {

class FileStream {
  std::string path_{};
  mutable std::ifstream handle_{};
  std::size_t file_size_{};

 public:
  FileStream() = default;
  explicit FileStream(std::string path);

  [[nodiscard]] const std::string &path() const noexcept;
  [[nodiscard]] const std::string &url() const noexcept;
  [[nodiscard]] std::size_t size() const;

  void seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  [[nodiscard]] std::size_t tellg() const noexcept;
  [[nodiscard]] bool eof() const noexcept;

  void read(std::string &buffer, std::size_t count);
  void read(char *buffer, std::size_t count);

  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  [[nodiscard]] T read();
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void read(T &buffer);

  template <typename Tin, typename Tout = std::make_signed_t<Tin>,
            typename std::enable_if<std::is_integral<Tin>::value>::type * = nullptr>
  [[nodiscard]] Tout read_as_signed();
  template <typename Tin, typename Tout = std::make_unsigned_t<Tin>,
            typename std::enable_if<std::is_integral<Tin>::value>::type * = nullptr>
  [[nodiscard]] Tout read_as_unsigned();
  template <typename Tin, typename std::enable_if<std::is_arithmetic<Tin>::value>::type * = nullptr>
  [[nodiscard]] double read_as_double();

  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  void read(std::vector<T> &buffer);
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
  [[nodiscard]] std::vector<T> read(std::size_t size);

  void append(std::string &buffer, std::size_t count);

  bool getline(std::string &buffer, char delim = '\n');
  [[nodiscard]] std::string getline(char delim = '\n');

 private:
  [[nodiscard]] std::streampos new_pos(std::streamoff offset, std::ios::seekdir way);
  [[nodiscard]] static std::ifstream open_file(const std::string &path,
                                               std::ifstream::openmode mode);
};
}  // namespace hictk::internal::filestream

#include "../../../filestream_impl.hpp"
