// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
#pragma comment(lib, "Rpcrt4.Lib")
#else
#error "Unsupported OS"
#endif

#include <atomic>
#include <cstdio>
#include <filesystem>
#include <utility>

namespace hictk::test {
// The point of this class is to provide a reliable way to create a directory that automatically
// deletes istelf and its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class SelfDeletingFolder {
  std::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};

 public:
  [[maybe_unused]] SelfDeletingFolder() {
    try {
      _path = create_uniq_temp_dir(std::filesystem::temp_directory_path());
    } catch (const std::filesystem::filesystem_error& e) {
      // Workaround spurious CI failures due to missing /tmp folder exception
      _path = create_uniq_temp_dir("test/data/unit_tests/scratch");
    }
  }

  [[maybe_unused]] explicit SelfDeletingFolder(std::filesystem::path path,
                                               bool delete_on_destruction = true)
      : _path(std::move(path)), _delete_on_destruction(delete_on_destruction) {
    std::filesystem::create_directories(_path);
  }

  [[maybe_unused]] explicit SelfDeletingFolder(bool delete_on_destruction) : SelfDeletingFolder() {
    this->set_delete_on_destruction(delete_on_destruction);
  }

  SelfDeletingFolder(const SelfDeletingFolder& other) = delete;
  SelfDeletingFolder(SelfDeletingFolder&& other) = delete;

  ~SelfDeletingFolder() {
    if (this->get_delete_on_destruction()) {
      std::filesystem::remove_all(this->_path);
    }
  }

  [[nodiscard]] const std::filesystem::path& operator()() const noexcept { return this->_path; }
  [[maybe_unused]] [[nodiscard]] bool get_delete_on_destruction() const noexcept {
    return this->_delete_on_destruction;
  }

  [[maybe_unused]] void set_delete_on_destruction(const bool flag) noexcept {
    this->_delete_on_destruction = flag;
  }

  SelfDeletingFolder& operator=(const SelfDeletingFolder& other) = delete;
  SelfDeletingFolder& operator=(SelfDeletingFolder&& other) = delete;

  [[nodiscard]] static std::filesystem::path create_uniq_temp_dir(
      const std::filesystem::path& tmpdir) {
#ifdef _WIN32
    UUID uid;
    std::filesystem::path dir{};

    do {
      auto status = UuidCreate(&uid);
      if (status != RPC_S_OK) {
        continue;
      }
      char* str;
      auto res = UuidToStringA(&uid, reinterpret_cast<RPC_CSTR*>(&str));
      if (res != RPC_S_OK) {
        continue;
      }

      dir = tmpdir / std::filesystem::path(std::string{"hictk-tmp-"} + std::string{str});
    } while (!std::filesystem::create_directories(dir));

    return dir;
#else
    return {mkdtemp((tmpdir / "hictk-tmp-XXXXXXXXXX").string().data())};
#endif
  }
};
}  // namespace hictk::test
