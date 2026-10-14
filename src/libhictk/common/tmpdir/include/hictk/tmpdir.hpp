// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>
#include <filesystem>

namespace hictk::internal {

// The point of this class is to provide a reliable way to create a directory that automatically
// deletes itself and its content by default.
// This can be prevented by setting the internal flag to true.
// The default constructor will create a unique, randomly named directory under the system tmpdir
class TmpDir {
  std::filesystem::path _path;
  std::atomic<bool> _delete_on_destruction{true};

  void delete_at_exit();

 public:
  [[maybe_unused]] TmpDir();
  [[maybe_unused]] explicit TmpDir(std::filesystem::path path);
  [[maybe_unused]] explicit TmpDir(const std::filesystem::path& prefix, bool delete_on_destruction);
  [[maybe_unused]] explicit TmpDir(bool delete_on_destruction);

  TmpDir(const TmpDir& other) = delete;
  TmpDir(TmpDir&& other) = delete;

  TmpDir& operator=(const TmpDir& other) = delete;
  TmpDir& operator=(TmpDir&& other) = delete;

  ~TmpDir() noexcept;

  [[nodiscard]] const std::filesystem::path& operator()() const noexcept;
  [[maybe_unused]] [[nodiscard]] bool get_delete_on_destruction() const noexcept;
  [[maybe_unused]] void set_delete_on_destruction(bool flag) noexcept;

  [[nodiscard]] static std::filesystem::path default_temp_directory_path();

  [[nodiscard]] static std::filesystem::path create_uniq_temp_dir(
      const std::filesystem::path& tmpdir);
};

}  // namespace hictk::internal
