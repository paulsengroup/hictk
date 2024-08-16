// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#elif defined(_WIN32)
#include <cstdio>
#include <random>
#include <string>
#endif

#include <fmt/format.h>
#include <fmt/std.h>

#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <utility>

namespace hictk::internal {
// The point of this class is to provide a reliable way to create a directory that automatically
// deletes istelf and its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class TmpDir {
  std::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};

  void delete_at_exit() {
    if (_delete_on_destruction) {
      std::filesystem::remove_all(_path);
    }
  }

 public:
  [[maybe_unused]] TmpDir() {
    try {
      _path = create_uniq_temp_dir(default_temp_directory_path());
    } catch (const std::filesystem::filesystem_error&) {
      const auto called_from_ci = std::getenv("HICTK_CI") != nullptr;  // NOLINT(*-mt-unsafe)
      if (!called_from_ci) {
        throw;
      }
      // Workaround spurious CI failures due to missing /tmp folder exception
      _path = create_uniq_temp_dir(std::filesystem::current_path());
    }
  }

  [[maybe_unused]] explicit TmpDir(std::filesystem::path path)
      : _path(std::move(path)), _delete_on_destruction(true) {
    if (std::filesystem::exists(_path)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to use path \"{}\" as TmpDir: folder already exists"),
                      _path.string()));
    }
    std::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit TmpDir(const std::filesystem::path& prefix, bool delete_on_destruction)
      : _path(create_uniq_temp_dir(prefix)), _delete_on_destruction(delete_on_destruction) {
    std::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit TmpDir(bool delete_on_destruction) : TmpDir() {
    set_delete_on_destruction(delete_on_destruction);
  }

  TmpDir(const TmpDir& other) = delete;
  TmpDir(TmpDir&& other) = delete;

  ~TmpDir() {
    if (get_delete_on_destruction()) {
      std::filesystem::remove_all(_path);
    }
  }

  [[nodiscard]] const std::filesystem::path& operator()() const noexcept { return _path; }
  [[maybe_unused]] [[nodiscard]] bool get_delete_on_destruction() const noexcept {
    return _delete_on_destruction;
  }

  [[maybe_unused]] void set_delete_on_destruction(const bool flag) noexcept {
    _delete_on_destruction = flag;
  }

  TmpDir& operator=(const TmpDir& other) = delete;
  TmpDir& operator=(TmpDir&& other) = delete;

  [[nodiscard]] static std::filesystem::path default_temp_directory_path() {
    try {
      return std::filesystem::temp_directory_path();
    } catch (const std::filesystem::filesystem_error& e) {
      if (e.path1().empty()) {
        throw std::filesystem::filesystem_error(
            "unable to safely determine the path where to store temporary files: please make sure "
            "the environment variable TMPDIR is defined and pointing to an existing folder",
            e.code());
      }
      throw std::filesystem::filesystem_error(
          fmt::format(
              FMT_STRING("unable to safely determine the path where to store temporary "
                         "files: temporary folder is set to \"{}\" but folder does not exists"),
              e.path1()),
          e.code());
    }
  }

  [[nodiscard]] static std::filesystem::path create_uniq_temp_dir(
      const std::filesystem::path& tmpdir) {
    if (!std::filesystem::exists(tmpdir)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to use path {} as TmpDir: path does not exists"), tmpdir));
    }
#ifdef _WIN32
    std::random_device rd;
    std::mt19937_64 rand_eng(rd());
    std::filesystem::path dir{};

    const auto lb = static_cast<int>('A');
    const auto ub = static_cast<int>('Z');
    do {
      std::string str{10, '\0'};

      for (std::size_t i = 0; i < str.size(); ++i) {
        str[i] = static_cast<char>(std::uniform_int_distribution<int>{lb, ub}(rand_eng));
      }

      dir = tmpdir / std::filesystem::path(std::string{"hictk-tmp-"} + std::string{str});
    } while (!std::filesystem::create_directories(dir));

    return dir;
#else
    auto dir = (tmpdir / "hictk-tmp-XXXXXXXXXX").string();
    if (!mkdtemp(dir.data())) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("unable to use path {} as TmpDir: failed to create a temporary folder"),
          tmpdir));
    }
    return {dir};
#endif
  }
};
}  // namespace hictk::internal
