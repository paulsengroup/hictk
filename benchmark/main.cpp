#include <catch2/catch_session.hpp>

#include "hictk/benchmark/benchmark_installers.hpp"

using namespace hictk::benchmark;

int main(int argc, char* argv[]) {
  register_cooler_cis_queries_benchmarks();
  register_cooler_gw_queries_benchmarks();
  register_cooler_trans_queries_benchmarks();

  register_file_cis_queries_benchmarks();
  register_file_gw_queries_benchmarks();
  register_file_trans_queries_benchmarks();

  register_hic_cis_queries_benchmarks();
  register_hic_gw_queries_benchmarks();
  register_hic_trans_queries_benchmarks();

  return Catch::Session().run(argc, argv);
}
