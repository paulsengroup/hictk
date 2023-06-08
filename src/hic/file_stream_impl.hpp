// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace hicxx::internal::filestream {

inline auto FileStream::forward_or_guess_stream_type(const std::string &url,
                                                     StreamType known_type) noexcept -> StreamType {
    if (known_type != StreamType::AUTO) {
        return known_type;
    }

    if (!!std::ifstream(url)) {
        return StreamType::LOCAL;
    }
    return StreamType::REMOTE;
}

inline FileStream::FileStream(std::string url, StreamType type, std::size_t chunk_size,
                              std::string agent)
    : FileStream(forward_or_guess_stream_type(url, type) == StreamType::LOCAL
                     ? FileStream::Local(std::move(url))
                     : FileStream::Remote(std::move(url), chunk_size, std::move(agent))) {}

inline FileStream::FileStream(internal::LocalFileStream &&fs) noexcept : stream_(std::move(fs)) {}
inline FileStream::FileStream(internal::RemoteFileStream &&fs) noexcept : stream_(std::move(fs)) {}

inline FileStream FileStream::Local(std::string path) {
    return FileStream{internal::LocalFileStream(std::move(path))};
}

inline FileStream FileStream::Remote(std::string path, std::size_t chunk_size, std::string agent) {
    return FileStream{internal::RemoteFileStream(std::move(path), chunk_size, std::move(agent))};
}

inline bool FileStream::is_local() const noexcept {
    return std::holds_alternative<internal::LocalFileStream>(this->stream_);
}

inline bool FileStream::is_remote() const noexcept { return !this->is_local(); }

inline const std::string &FileStream::url() const noexcept {
    return std::visit([](const auto &fs) -> const std::string & { return fs.url(); },
                      this->stream_);
}

inline std::size_t FileStream::size() const noexcept {
    return std::visit([](const auto &fs) { return fs.size(); }, this->stream_);
}

inline void FileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
    return std::visit([&](auto &fs) { return fs.seekg(offset, way); }, this->stream_);
}

inline std::size_t FileStream::tellg() const noexcept {
    return std::visit([](const auto &fs) { return fs.tellg(); }, this->stream_);
}

inline bool FileStream::eof() const noexcept {
    return std::visit([](const auto &fs) { return fs.eof(); }, this->stream_);
}

inline void FileStream::read(std::string &buffer, std::size_t count) {
    return std::visit([&](auto &fs) { return fs.read(buffer, count); }, this->stream_);
}

inline void FileStream::read(char *buffer, std::size_t count) {
    return std::visit([&](auto &fs) { return fs.read(buffer, count); }, this->stream_);
}

inline void FileStream::append(std::string &buffer, std::size_t count) {
    return std::visit([&](auto &fs) { return fs.append(buffer, count); }, this->stream_);
}

inline bool FileStream::getline(std::string &buffer, char delim) {
    return std::visit([&](auto &fs) { return fs.getline(buffer, delim); }, this->stream_);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::read(T &buffer) {
    return std::visit([&](auto &fs) { return fs.read(buffer); }, this->stream_);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T FileStream::read() {
    T buffer{};
    read(buffer);
    return buffer;
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::read(std::vector<T> &buffer) {
    return std::visit([&](auto &fs) { return fs.read(buffer); }, this->stream_);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline std::vector<T> FileStream::read(std::size_t size) {
    return std::visit([&](auto &fs) { return fs.template read<T>(size); }, this->stream_);
}

inline std::string FileStream::getline(char delim) {
    return std::visit([&](auto &fs) { return fs.getline(delim); }, this->stream_);
}

}  // namespace hicxx::internal::filestream
