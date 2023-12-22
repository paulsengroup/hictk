// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cstdio>

namespace std {
template <>
struct default_delete<FILE> {
  void operator()(FILE* file) const { std::fclose(file); }  // NOLINT
};
}  // namespace std
