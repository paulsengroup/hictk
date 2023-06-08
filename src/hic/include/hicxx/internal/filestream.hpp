// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICXX_USE_CURL
#include <curl/curl.h>
#endif

#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include "hicxx/internal/common.hpp"

namespace hicxx::internal::filestream {

namespace internal {
#ifdef HICXX_USE_CURL
class RemoteFileStream {
   public:
    using CURL_ptr = UniquePtrWithDeleter<CURL>;

   private:
    std::string url_{};
    CURL_ptr handle_{};
    std::string buffer_{};
    std::size_t chunk_offset_{};
    std::size_t stream_pos_{};
    std::size_t stream_size_{};

   public:
    RemoteFileStream() = default;
    explicit RemoteFileStream(std::string url, std::size_t chunk_size = 64 * 1024,
                              std::string agent = "hicxx");

    [[nodiscard]] const std::string &url() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

    void seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
    [[nodiscard]] std::size_t tellg() const noexcept;
    [[nodiscard]] bool eof() const noexcept;

    void read(std::string &buffer, std::size_t count);
    void read(char *buffer, std::size_t count);

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    [[nodiscard]] T read();
    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    void read(T &buffer);

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    void read(std::vector<T> &buffer);
    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    [[nodiscard]] std::vector<T> read(std::size_t size);

    void append(std::string &buffer, std::size_t count);

    bool getline(std::string &buffer, char delim = '\n');
    [[nodiscard]] std::string getline(char delim = '\n');

   private:
    [[nodiscard]] std::streampos new_pos(std::streamoff offset, std::ios::seekdir way);
    [[nodiscard]] std::size_t available_bytes() const noexcept;

    [[nodiscard]] std::size_t first_chunk_pos() const noexcept;
    [[nodiscard]] std::size_t last_chunk_pos() const noexcept;

    [[nodiscard]] char *chunk_begin() noexcept;
    [[nodiscard]] const char *chunk_begin() const noexcept;

    [[nodiscard]] char *chunk_current() noexcept;
    [[nodiscard]] const char *chunk_current() const noexcept;

    [[nodiscard]] char *chunk_end() noexcept;
    [[nodiscard]] const char *chunk_end() const noexcept;

    void fetch_next_chunk();
    [[nodiscard]] std::size_t eof_pos() const noexcept;
    void mark_eof() noexcept;
    // void reserve_buffer(std::size_t new_size);
    static std::size_t write_memory_callback(void *contents, std::size_t size, std::size_t nmemb,
                                             void *userp) noexcept;
    [[nodiscard]] static auto init_CURL(const std::string &url, const std::string &agent)
        -> CURL_ptr;
    [[nodiscard]] static std::size_t get_stream_size(const std::string &url,
                                                     const std::string &agent);
};
#endif

class LocalFileStream {
    std::string path_{};
    mutable std::ifstream handle_{};
    std::size_t file_size_{};

   public:
    LocalFileStream() = default;
    explicit LocalFileStream(std::string path);

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
}  // namespace internal

#ifndef HICXX_USE_CURL
using FileStream = internal::LocalFileStream;
#else

class FileStream {
    enum class StreamType { AUTO, LOCAL, REMOTE };
    using StreamVariant = std::variant<internal::LocalFileStream, internal::RemoteFileStream>;

    StreamVariant stream_{internal::LocalFileStream()};

   public:
    FileStream() = default;
    explicit FileStream(std::string url, StreamType type = StreamType::AUTO,
                        std::size_t chunk_size = 64 * 1024, std::string agent = "filestream");
    static FileStream Local(std::string path);
    static FileStream Remote(std::string path, std::size_t chunk_size = 64 * 1024,
                             std::string agent = "hicxx");

    bool is_local() const noexcept;
    bool is_remote() const noexcept;

    const std::string &url() const noexcept;
    std::size_t size() const noexcept;

    void seekg(std::streamoff offset, std::ios::seekdir way = std::ios::beg);
    std::size_t tellg() const noexcept;
    bool eof() const noexcept;

    void read(std::string &buffer, std::size_t count);
    void read(char *buffer, std::size_t count);

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    T read();
    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    void read(T &buffer);

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    void read(std::vector<T> &buffer);
    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    std::vector<T> read(std::size_t size);

    void append(std::string &buffer, std::size_t count);

    bool getline(std::string &buffer, char delim = '\n');
    std::string getline(char delim = '\n');

   private:
    static StreamType forward_or_guess_stream_type(const std::string &url,
                                                   StreamType known_type) noexcept;
    explicit FileStream(internal::LocalFileStream &&fs) noexcept;
    explicit FileStream(internal::RemoteFileStream &&fs) noexcept;
};
#endif

}  // namespace hicxx::internal::filestream

#include "../../../local_fstream_impl.hpp"
#ifdef HICXX_USE_CURL
#include "../../../file_stream_impl.hpp"
#include "../../../remote_fstream_impl.hpp"
#endif
