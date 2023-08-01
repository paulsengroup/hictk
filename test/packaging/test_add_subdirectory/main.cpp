// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <iostream>

#include "hictk/version.hpp"

int main() { std::cout << hictk::config::version::str_long() << "\n"; }
