// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace hictk::filestream {

template <typename Mutex = std::mutex>
class FileStream {
  static_assert(sizeof(char) == 1, "char must be 1 byte wide!");
  static_assert(sizeof(float) == 4, "float must be 4 bytes wide!");
  static_assert(sizeof(double) == 8, "double must be 8 bytes wide!");
  std::string _path{};
  mutable std::shared_ptr<Mutex> _mtx{};
  mutable std::ifstream _ifs{};
  mutable std::ofstream _ofs{};
  std::streamsize _file_size{};

 public:
  FileStream() = default;
  // pass nullptr as mtx to disable file locking
  FileStream(std::string path, std::shared_ptr<Mutex> mtx, std::ios::openmode mode = std::ios::in);
  static FileStream create(std::string path, std::shared_ptr<Mutex> mtx);

  FileStream(const FileStream &other) = delete;
  FileStream(FileStream &&other) noexcept = default;

  ~FileStream() noexcept = default;

  FileStream &operator=(const FileStream &other) = delete;
  FileStream &operator=(FileStream &&other) noexcept = default;

  [[nodiscard]] const std::string &path() const noexcept;

  // seek*
  void seekg(std::streampos position);
  void unsafe_seekg(std::streampos position);
  void seekp(std::streampos position);
  void unsafe_seekp(std::streampos position);

  void seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  void unsafe_seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  void seekp(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
  void unsafe_seekp(std::streamoff offset, std::ios::seekdir way = std::ios::beg);

  // tell*
  [[nodiscard]] std::streampos tellg() const;
  [[nodiscard]] std::streampos unsafe_tellg() const;
  [[nodiscard]] std::streampos tellp() const;
  [[nodiscard]] std::streampos unsafe_tellp() const;

  // others
  [[nodiscard]] std::streamsize size() const;
  [[nodiscard]] std::streamsize unsafe_size() const;

  [[nodiscard]] bool eof() const;
  [[nodiscard]] bool unsafe_eof() const;

  void flush();
  void unsafe_flush();

  [[nodiscard]] static std::string get_underlying_os_error();

  // Locks mutex protecting the underlying file streams.
  // IMPORTANT: while the lock is held, only unsafe_* methods can be used.
  [[nodiscard]] std::unique_lock<Mutex> lock() const;
  [[nodiscard]] bool is_locked() const noexcept;

  // read char*
  void read(char *buffer, std::size_t count);
  // this method (as well as all other seek_* methods) return a pair with the offset before seek and
  // the offset after read/write
  std::pair<std::streampos, std::streampos> seek_and_read(std::streamoff offset, char *buffer,
                                                          std::size_t count,
                                                          std::ios::seekdir way = std::ios::beg);
  void unsafe_read(char *buffer, std::size_t count);

  // read std::string
  [[nodiscard]] std::string read(std::size_t count);
  void read(std::string &buffer, std::size_t count);
  void read_append(std::string &buffer, std::size_t count);
  std::pair<std::streampos, std::streampos> seek_and_read(std::streamoff offset,
                                                          std::string &buffer, std::size_t count,
                                                          std::ios::seekdir way = std::ios::beg);

  // getline
  bool getline(std::string &buffer, char delim = '\n');
  [[nodiscard]] std::string getline(char delim = '\n');
  std::tuple<bool, std::streampos, std::streampos> seek_and_getline(
      std::streamoff offset, std::string &buffer, std::ios::seekdir way = std::ios::beg,
      char delim = '\n');
  bool unsafe_getline(std::string &buffer, char delim);

  // read T
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  [[nodiscard]] T read();
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void read(T &buffer);
  template <typename Tin, typename Tout = std::make_signed_t<Tin>,
            typename std::enable_if_t<std::is_integral_v<Tin>> * = nullptr>
  [[nodiscard]] Tout read_as_signed();
  template <typename Tin, typename Tout = std::make_unsigned_t<Tin>,
            typename std::enable_if_t<std::is_integral_v<Tin>> * = nullptr>
  [[nodiscard]] Tout read_as_unsigned();
  template <typename Tin, typename std::enable_if_t<std::is_arithmetic_v<Tin>> * = nullptr>
  [[nodiscard]] double read_as_double();
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void unsafe_read(T &buffer);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  [[nodiscard]] T unsafe_read();

  // read vector<T>
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void read(std::vector<T> &buffer);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  [[nodiscard]] std::vector<T> read_vector(std::size_t size);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  std::pair<std::streampos, std::streampos> seek_and_read(std::streamoff offset,
                                                          std::vector<T> &buffer,
                                                          std::ios::seekdir way = std::ios::beg);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void unsafe_read(std::vector<T> &buffer);

  // write const char*
  void write(const char *buffer, std::size_t count);
  std::pair<std::streampos, std::streampos> seek_and_write(std::streamoff offset,
                                                           const char *buffer, std::size_t count,
                                                           std::ios::seekdir way = std::ios::beg);
  void unsafe_write(const char *buffer, std::size_t count);
  std::pair<std::streampos, std::streampos> append(const char *buffer, std::size_t count);

  // write std::string
  void write(std::string_view buffer);
  std::pair<std::streampos, std::streampos> seek_and_write(std::streamoff offset,
                                                           std::string_view buffer,
                                                           std::ios::seekdir way = std::ios::beg);
  std::pair<std::streampos, std::streampos> append(std::string_view buffer);

  // write T
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void write(T buffer);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void unsafe_write(T buffer);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  std::pair<std::streampos, std::streampos> seek_and_write(std::streamoff offset, T buffer,
                                                           std::ios::seekdir way = std::ios::beg);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  std::pair<std::streampos, std::streampos> append(T buffer);

  // write vector<T>
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  void write(const std::vector<T> &buffer);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  std::pair<std::streampos, std::streampos> seek_and_write(std::streamoff offset,
                                                           const std::vector<T> &buffer,
                                                           std::ios::seekdir way = std::ios::beg);
  template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> * = nullptr>
  std::pair<std::streampos, std::streampos> append(const std::vector<T> &buffer);

 private:
  [[nodiscard]] std::streampos new_posg(std::streamoff offset, std::ios::seekdir way);
  [[nodiscard]] std::streampos new_posg_checked(std::streamoff offset, std::ios::seekdir way);
  [[nodiscard]] std::streampos new_posp(std::streamoff offset, std::ios::seekdir way);
  void update_file_size();
  void unsafe_update_file_size();
  [[nodiscard]] static std::ifstream open_file_read(const std::string &path,
                                                    std::ifstream::openmode mode);
  [[nodiscard]] static std::ofstream open_file_write(const std::string &path,
                                                     std::ofstream::openmode mode);
  [[nodiscard]] static std::string get_underlying_os_error(int errno_);
  static void get_underlying_os_error(int errno_, std::string &buffer);

  template <typename T, typename I1, typename I2>
  static void validate_read(I1 bytes_read, I2 count_expected,
                            std::string_view prefix = "FileStream::read*()");
};
}  // namespace hictk::filestream

#include "./impl/filestream_impl.hpp"  // NOLINT
