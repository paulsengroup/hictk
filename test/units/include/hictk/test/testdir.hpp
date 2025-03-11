// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

#include "hictk/tmpdir.hpp"

namespace hictk::test {

inline const internal::TmpDir testdir{true};              // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)

}  // namespace hictk::test
