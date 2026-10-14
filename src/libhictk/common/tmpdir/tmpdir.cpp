// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tmpdir.hpp"

#ifdef HICTK_MKDTEMP_AVAIL
#include <unistd.h>
#endif

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <atomic>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace hictk::internal {

void TmpDir::delete_at_exit() {
  if (_delete_on_destruction) {
    std::filesystem::remove_all(_path);  // NOLINT
  }
}

[[nodiscard]] static std::filesystem::path ctor_helper() {
  try {
    return TmpDir::create_uniq_temp_dir(TmpDir::default_temp_directory_path());
  } catch (const std::filesystem::filesystem_error&) {
    const auto called_from_ci = std::getenv("HICTK_CI") != nullptr;  // NOLINT(*-mt-unsafe)
    if (!called_from_ci) {
      throw;
    }
    // Workaround spurious CI failures due to missing /tmp folder exception
    return TmpDir::create_uniq_temp_dir(std::filesystem::current_path());
  }
}

TmpDir::TmpDir() : _path(ctor_helper()), _delete_on_destruction(true) {}

TmpDir::TmpDir(std::filesystem::path path) : _path(std::move(path)), _delete_on_destruction(true) {
  if (std::filesystem::exists(_path)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to use path \"{}\" as TmpDir: folder already exists"), _path.string()));
  }
  std::filesystem::create_directories(_path);
}

TmpDir::TmpDir(const std::filesystem::path& prefix, bool delete_on_destruction)
    : _path(create_uniq_temp_dir(prefix)), _delete_on_destruction(delete_on_destruction) {
  std::filesystem::create_directories(_path);
}

TmpDir::TmpDir(bool delete_on_destruction) : TmpDir() {
  set_delete_on_destruction(delete_on_destruction);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
TmpDir::~TmpDir() noexcept {
  try {
    delete_at_exit();
  } catch (std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to delete temporary folder \"{}\": {}"), _path, e.what());
  } catch (...) {
    SPDLOG_WARN(FMT_STRING("failed to delete temporary folder \"{}\""), _path);
  }
}

const std::filesystem::path& TmpDir::operator()() const noexcept { return _path; }
bool TmpDir::get_delete_on_destruction() const noexcept { return _delete_on_destruction.load(); }

void TmpDir::set_delete_on_destruction(const bool flag) noexcept { _delete_on_destruction = flag; }

std::filesystem::path TmpDir::default_temp_directory_path() {
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
                       "files: temporary folder is set to \"{}\" but folder does not exist"),
            e.path1()),
        e.code());
  }
}

[[nodiscard]] static std::filesystem::path create_uniq_temp_dir_fallback(
    const std::filesystem::path& tmpdir, std::size_t max_attempts = 100) {
  assert(max_attempts > 0);

  std::random_device rd;
  std::mt19937_64 rand_eng(rd());
  std::filesystem::path dir{};

  const auto lb = static_cast<int>('A');
  const auto ub = static_cast<int>('Z');
  std::string str{10, '\0'};
  for (std::uint8_t i = 0; i < max_attempts; ++i) {
    std::generate(str.begin(), str.end(), [&]() {
      return static_cast<char>(std::uniform_int_distribution<int>{lb, ub}(rand_eng));
    });

    dir = tmpdir / std::filesystem::path(fmt::format(FMT_STRING("hictk-tmp-{}"), str));
    if (std::filesystem::create_directories(dir)) {
      return dir;
    }
  }

  throw std::runtime_error(fmt::format(
      FMT_STRING("unable to use path {} as TmpDir: failed to create a temporary folder"), tmpdir));
}

[[nodiscard]] static std::filesystem::path create_uniq_temp_dir_mkdtemp(
    const std::filesystem::path& tmpdir) {
  auto dir = (tmpdir / "hictk-tmp-XXXXXXXXXX").string();
  if (mkdtemp(dir.data())) {
    return {dir};
  }
  return {};
}

std::filesystem::path TmpDir::create_uniq_temp_dir(const std::filesystem::path& tmpdir) {
  if (!std::filesystem::exists(tmpdir)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to use path {} as TmpDir: path does not exists"), tmpdir));
  }

#ifdef HICTK_MKDTEMP_AVAIL
  if (auto dir = create_uniq_temp_dir_mkdtemp(tmpdir); !dir.empty()) {
    return dir;
  }
  SPDLOG_DEBUG(
      "TmpDir::create_uniq_temp_dir: failed to create temporary directory using mkdtemp: using "
      "fallback method...");
#endif

  return create_uniq_temp_dir_fallback(tmpdir);
}
}  // namespace hictk::internal
