// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

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
