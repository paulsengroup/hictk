// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictk/tmpdir.hpp"

namespace hictk {

namespace test {
inline const internal::TmpDir testdir{true};  // NOLINT(cert-err58-cpp)

}  // namespace test

namespace cooler::test::attribute {

static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::attribute

namespace cooler::test::balancing {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::balancing

namespace cooler::test::cooler_file {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::cooler_file

namespace cooler::test::multires_cooler_file {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::multires_cooler_file

namespace cooler::test::singlecell_cooler_file {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::singlecell_cooler_file

namespace cooler::test::dataset {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::dataset

namespace cooler::test::pixel_selector {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::pixel_selector

namespace filestream::test {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace filestream::test

namespace hic::test::file_reader {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::file_reader

namespace hic::test::file_writer {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::file_writer

namespace hic::test::utils {
static const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::utils

}  // namespace hictk
