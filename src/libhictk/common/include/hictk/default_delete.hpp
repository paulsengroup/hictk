// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <libdeflate.h>
#include <zstd.h>

#include <memory>

template <>
struct std::default_delete<ZSTD_CCtx_s> {
  void operator()(ZSTD_CCtx_s* ctx) const { ZSTD_freeCCtx(ctx); }  // NOLINT
};

template <>
struct std::default_delete<ZSTD_DCtx_s> {
  void operator()(ZSTD_DCtx_s* ctx) const { ZSTD_freeDCtx(ctx); }  // NOLINT
};

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
};

template <>
struct std::default_delete<FILE> {
  void operator()(FILE* file) const { std::fclose(file); }  // NOLINT
};
