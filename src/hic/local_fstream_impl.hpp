// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hicxx::internal::filestream::internal {

inline LocalFileStream::LocalFileStream(std::string path)
    : path_(std::move(path)),
      handle_(open_file(path_, std::ios::binary | std::ios::ate)),
      file_size_(static_cast<std::size_t>(handle_.tellg())) {
    handle_.seekg(std::ios::beg);
}

inline const std::string &LocalFileStream::path() const noexcept { return this->path_; }
inline const std::string &LocalFileStream::url() const noexcept { return this->path(); }

inline std::size_t LocalFileStream::size() const { return this->file_size_; }

inline void LocalFileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
    const auto new_pos = this->new_pos(offset, way);
    if (new_pos < 0 || new_pos >= std::int64_t(this->size() + 1)) {
        throw std::runtime_error("caught an attempt of out-of-bound read");
    }
    this->handle_.seekg(new_pos, std::ios::beg);
}

inline std::size_t LocalFileStream::tellg() const noexcept {
    return static_cast<std::size_t>(this->handle_.tellg());
}

inline bool LocalFileStream::eof() const noexcept { return this->handle_.eof(); }

inline void LocalFileStream::read(std::string &buffer, std::size_t count) {
    buffer.resize(count);
    if (count > 0) {
        return this->read(&buffer.front(), count);
    }
}

inline void LocalFileStream::read(char *buffer, std::size_t count) {
    this->handle_.read(buffer, std::int64_t(count));
}

inline void LocalFileStream::append(std::string &buffer, std::size_t count) {
    const auto buff_size = buffer.size();
    buffer.resize(buffer.size() + count);

    this->handle_.read(&(*buffer.begin()) + buff_size, std::int64_t(count));
}

inline bool LocalFileStream::getline(std::string &buffer, char delim) {
    buffer.clear();
    if (this->eof()) {
        this->handle_.setstate(std::ios::badbit);
    }
    try {
        return !!std::getline(this->handle_, buffer, delim);
    } catch (const std::exception &e) {
        if (this->handle_.eof() && !this->handle_.bad()) {
            return !!this->handle_;
        }
        throw;
    }
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void LocalFileStream::read(T &buffer) {
    static_assert(sizeof(char) == 1, "");
    return this->read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T LocalFileStream::read() {
    T buffer{};
    this->read(buffer);
    return buffer;
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void LocalFileStream::read(std::vector<T> &buffer) {
    static_assert(sizeof(char) == 1, "");
    return this->read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline std::vector<T> LocalFileStream::read(std::size_t size) {
    assert(size != 0);
    std::vector<T> buffer(size);
    this->read(buffer);
    return buffer;
}

inline std::string LocalFileStream::getline(char delim) {
    std::string buffer{};
    this->getline(buffer, delim);
    return buffer;
}

inline std::streampos LocalFileStream::new_pos(std::streamoff offset, std::ios::seekdir way) {
    switch (way) {
        case std::ios::beg:
            return static_cast<std::streampos>(offset);
        case std::ios::cur:
            return std::int64_t(this->tellg()) + offset;
        case std::ios::end:
            return std::int64_t(this->file_size_) + offset;
        default:
            assert(false);
            std::abort();
    }
}

inline std::ifstream LocalFileStream::open_file(const std::string &path,
                                                std::ifstream::openmode mode) {
    std::ifstream ifs;
    ifs.exceptions(ifs.exceptions() | std::ios::failbit | std::ios::badbit);
    ifs.open(path, mode);
    return ifs;
}

}  // namespace hicxx::internal::filestream::internal
