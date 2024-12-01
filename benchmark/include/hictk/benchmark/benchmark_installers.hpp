// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::benchmark {

void register_cooler_cis_queries_benchmarks();
void register_cooler_gw_queries_benchmarks();
void register_cooler_trans_queries_benchmarks();

void register_file_cis_queries_benchmarks();
void register_file_gw_queries_benchmarks();
void register_file_trans_queries_benchmarks();

void register_hic_cis_queries_benchmarks();
void register_hic_gw_queries_benchmarks();
void register_hic_trans_queries_benchmarks();

}  // namespace hictk::benchmark
