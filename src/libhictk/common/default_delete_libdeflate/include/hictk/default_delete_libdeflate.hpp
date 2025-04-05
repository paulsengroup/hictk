// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <libdeflate.h>

#include <memory>

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
};
