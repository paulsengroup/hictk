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

inline RemoteFileStream::RemoteFileStream(std::string url, std::size_t chunk_size,
                                          std::string agent)
    : url_(std::move(url)),
      handle_(init_CURL(url_, agent)),
      buffer_(chunk_size, '\0'),
      stream_size_(get_stream_size(url_, agent)) {
    buffer_.clear();
}

inline std::size_t RemoteFileStream::write_memory_callback(void *contents, std::size_t size,
                                                           std::size_t nmemb,
                                                           void *userp) noexcept {
    assert(userp);
    assert(contents);
    std::size_t realsize = size * nmemb;
    auto &buffer = *reinterpret_cast<std::string *>(userp);  // NOLINT
    try {
        buffer.append(static_cast<const char *>(contents), realsize);
    } catch (const std::bad_alloc &e) {
        return 0;
    }
    return realsize;
}

inline auto RemoteFileStream::init_CURL(const std::string &url, const std::string &agent)
    -> CURL_ptr {
    auto curl = CURL_ptr(curl_easy_init(), &curl_easy_cleanup);
    if (curl.get()) {
        curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION,
                         RemoteFileStream::write_memory_callback);
        curl_easy_setopt(curl.get(), CURLOPT_URL, &url.front());
        curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, agent.c_str());
    } else {
        throw std::runtime_error("Unable to initialize curl");
    }
    return curl;
}

inline std::size_t RemoteFileStream::get_stream_size(const std::string &url,
                                                     const std::string &agent) {
    auto discardData = +[]([[maybe_unused]] void *buffer, std::size_t size, std::size_t nmemb,
                           [[maybe_unused]] void *userp) -> std::size_t { return size * nmemb; };

    auto curl = CURL_ptr(curl_easy_init(), &curl_easy_cleanup);
    if (!curl.get()) {
        throw std::runtime_error("Unable to initialize curl");
    }
    curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl.get(), CURLOPT_HEADER, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_NOBODY, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, discardData);
    curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, agent.c_str());

    auto res = curl_easy_perform(curl.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("Unable to fetch metadata for " + url + ": " +
                                 std::string(curl_easy_strerror(res)));
    }
    curl_off_t cl{};
    res = curl_easy_getinfo(curl.get(), CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &cl);
    if (res != CURLE_OK || cl == -1) {
        const auto reason =
            res == CURLE_OK ? std::string{} : ": " + std::string(curl_easy_strerror(res));
        throw std::runtime_error("Unable to fetch content length for " + url + reason);
    }
    if (cl == -1) {
        throw std::runtime_error("Unable to fetch content length for " + url);
    }
    return static_cast<std::size_t>(cl);
}

inline const std::string &RemoteFileStream::url() const noexcept { return this->url_; }
inline std::size_t RemoteFileStream::size() const noexcept { return this->stream_size_; }

inline void RemoteFileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
    const auto new_pos = [&]() {
        const auto p = this->new_pos(offset, way);
        if (p < 0 || static_cast<std::size_t>(p) >= this->eof_pos()) {
            throw std::runtime_error("caught an attempt of out-of-bound read");
        }
        return static_cast<std::size_t>(p);
    }();

    if (new_pos < this->first_chunk_pos() || new_pos >= this->last_chunk_pos()) {
        this->buffer_.clear();
        this->chunk_offset_ = 0;
        this->stream_pos_ = new_pos;
    } else {
        this->chunk_offset_ = new_pos - this->stream_pos_;
    }
}

inline std::size_t RemoteFileStream::tellg() const noexcept {
    const auto pos = this->stream_pos_ + this->chunk_offset_;
    assert(pos <= this->eof_pos());
    return pos;
}

inline bool RemoteFileStream::eof() const noexcept { return this->tellg() == this->eof_pos(); }

inline void RemoteFileStream::read(std::string &buffer, std::size_t count) {
    buffer.resize(count);
    if (count > 0) {
        return this->read(&buffer.front(), count);
    }
}

inline void RemoteFileStream::read(char *buffer, std::size_t count) {
    if (this->tellg() + count > this->size()) {
        throw std::runtime_error("caught an attempt of out-of-bound read");
    }
    if (count == 0) {
        return;
    }

    if (this->available_bytes() == 0) {
        this->fetch_next_chunk();
    }

    if (count <= this->available_bytes()) {
        std::copy(this->chunk_current(), this->chunk_current() + count, buffer);
        this->chunk_offset_ += count;
        return;
    }

    std::copy(this->chunk_current(), this->chunk_end(), buffer);
    count -= this->available_bytes();
    buffer += this->available_bytes();
    this->chunk_offset_ = this->buffer_.size();
    this->read(buffer, count);
}

inline void RemoteFileStream::append(std::string &buffer, std::size_t count) {
    if (this->tellg() + count == this->size() + 1) {
        this->mark_eof();
        return;
    }
    if (this->tellg() + count > this->size()) {
        throw std::runtime_error("caught an attempt of out-of-bound read");
    }

    if (count == 0) {
        return;
    }

    if (this->available_bytes() == 0) {
        this->fetch_next_chunk();
    }

    if (count <= this->available_bytes()) {
        buffer.append(this->chunk_current(), count);
        this->chunk_offset_ += count;
        return;
    }

    buffer.append(this->chunk_current(), this->available_bytes());

    count -= this->available_bytes();
    this->chunk_offset_ = this->buffer_.size();
    this->append(buffer, count - this->available_bytes());
}

inline bool RemoteFileStream::getline(std::string &buffer, char delim) {
    buffer.clear();
    if (this->eof()) {
        throw std::runtime_error("caught an attempt of out-of-bound read");
    }

    while (!this->eof()) {
        const auto eol_pos = this->buffer_.find(delim, this->chunk_offset_);
        if (eol_pos != std::string::npos) {
            if (eol_pos != this->chunk_offset_) {
                const auto *first = this->chunk_current();
                const auto *last = this->chunk_begin() + eol_pos;
                buffer.append(first, last);
            }
            this->chunk_offset_ = eol_pos + 1;
            return !this->eof();
        }
        if (this->available_bytes() != 0) {
            buffer.append(this->chunk_current(), this->available_bytes());
        }
        this->chunk_offset_ = this->buffer_.size();

        this->fetch_next_chunk();
    }
    return !this->eof();
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void RemoteFileStream::read(T &buffer) {
    static_assert(sizeof(char) == 1, "");
    return this->read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T RemoteFileStream::read() {
    T buffer{};
    this->read(buffer);
    return buffer;
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void RemoteFileStream::read(std::vector<T> &buffer) {
    static_assert(sizeof(char) == 1, "");
    return this->read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline std::vector<T> RemoteFileStream::read(std::size_t size) {
    assert(size != 0);
    std::vector<T> buffer(size);
    this->read(buffer);
    return buffer;
}

inline std::string RemoteFileStream::getline(char delim) {
    std::string buffer{};
    this->getline(buffer, delim);
    return buffer;
}

inline std::streampos RemoteFileStream::new_pos(std::streamoff offset, std::ios::seekdir way) {
    switch (way) {
        case std::ios::beg:
            return static_cast<std::streampos>(offset);
        case std::ios::cur:
            return std::streampos(std::int64_t(this->tellg())) + offset;
        case std::ios::end:
            return std::streampos(std::int64_t(this->stream_size_)) + offset;
        default:
            assert(false);
            std::abort();
    }
}

inline std::size_t RemoteFileStream::available_bytes() const noexcept {
    assert(this->chunk_offset_ <= this->buffer_.size());
    return this->buffer_.size() - this->chunk_offset_;
}

inline std::size_t RemoteFileStream::first_chunk_pos() const noexcept { return this->stream_pos_; }

inline std::size_t RemoteFileStream::last_chunk_pos() const noexcept {
    return this->first_chunk_pos() + this->buffer_.size();
}

inline char *RemoteFileStream::chunk_begin() noexcept { return &(*this->buffer_.begin()); }

inline const char *RemoteFileStream::chunk_begin() const noexcept { return this->buffer_.data(); }

inline char *RemoteFileStream::chunk_current() noexcept {
    return this->chunk_begin() + this->chunk_offset_;
}

inline const char *RemoteFileStream::chunk_current() const noexcept {
    return this->chunk_begin() + this->chunk_offset_;
}

inline char *RemoteFileStream::chunk_end() noexcept {
    return this->chunk_begin() + this->buffer_.size();
}

inline const char *RemoteFileStream::chunk_end() const noexcept {
    return this->chunk_begin() + this->buffer_.size();
}

inline void RemoteFileStream::fetch_next_chunk() {
    if (this->eof()) {
        throw std::runtime_error("caught an attempt of out-of-bound read");
    }
    if (this->tellg() == this->size()) {  // EOF
        this->mark_eof();
        return;
    }
    assert(this->handle_);
    const auto first_pos = this->tellg();
    const auto last_pos = std::min(this->stream_size_, first_pos + this->buffer_.capacity() - 1);
    const auto range = std::to_string(first_pos) + "-" + std::to_string(last_pos);

    curl_easy_setopt(this->handle_.get(), CURLOPT_WRITEDATA,
                     reinterpret_cast<const void *>(&this->buffer_));
    curl_easy_setopt(this->handle_.get(), CURLOPT_RANGE, range.c_str());
    this->buffer_.clear();
    auto res = curl_easy_perform(this->handle_.get());
    if (res != CURLE_OK) {
        throw std::runtime_error("curl_easy_perform failed: " +
                                 std::string(curl_easy_strerror(res)));
    }

    this->chunk_offset_ = 0;
    this->stream_pos_ = first_pos;
}

inline std::size_t RemoteFileStream::eof_pos() const noexcept { return this->stream_size_ + 1; }

inline void RemoteFileStream::mark_eof() noexcept {
    this->buffer_.clear();
    this->chunk_offset_ = 0;
    this->stream_pos_ = this->eof_pos();
}
}  // namespace hicxx::internal::filestream::internal
