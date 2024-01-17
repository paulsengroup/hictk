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

inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::attribute

namespace cooler::test::balancing {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::balancing

namespace cooler::test::cooler_file {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::cooler_file

namespace cooler::test::multires_cooler_file {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::multires_cooler_file

namespace cooler::test::singlecell_cooler_file {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::singlecell_cooler_file

namespace cooler::test::dataset {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::dataset

namespace cooler::test::pixel_selector {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace cooler::test::pixel_selector

namespace hic::test::filestream {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::filestream

namespace hic::test::file_reader {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::file_reader

namespace hic::test::file_writer {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::file_writer

namespace hic::test::utils {
inline const auto& testdir = hictk::test::testdir;
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hic::test::file_writer

}  // namespace hictk
