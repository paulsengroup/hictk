// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/numpy.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstdint>
#include <variant>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::fuzzer {

template <typename N>
using NumpyArray = pybind11::array_t<N, pybind11::array::c_style | pybind11::array::forcecast>;

template <typename N>
using Eigen2DDense = Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename N>
using EigenSparse = Eigen::SparseMatrix<N, Eigen::RowMajor>;

// clang-format off
using PixelBuffer =
    std::variant<
        std::vector<ThinPixel<std::int32_t>>,
        std::vector<ThinPixel<double>>,
        std::vector<Pixel<std::int32_t>>,
        std::vector<Pixel<double>>>;
// clang-format on

}  // namespace hictk::fuzzer
