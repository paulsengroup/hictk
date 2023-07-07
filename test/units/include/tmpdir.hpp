// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};  // NOLINT(cert-err58-cpp)

}  // namespace hictk::test

namespace hictk::cooler::test::cooler_file {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::cooler::test::cooler_file
