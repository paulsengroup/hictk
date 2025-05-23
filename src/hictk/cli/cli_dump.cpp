// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "hictk/cooler/validation.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/validators.hpp"

namespace hictk::tools {

void Cli::make_dump_subcommand() {
  auto& sc = *_cli.add_subcommand("dump",
                                  "Read interactions and other kinds of data from .hic and Cooler "
                                  "files and write them to stdout.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(_config.index() == 0);
                    _config = DumpConfig{};
                  });

  _config = DumpConfig{};
  auto& c = std::get<DumpConfig>(_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic, .cool or .mcool file (Cooler URI syntax supported).")
      ->check(IsValidCoolerFile | IsValidHiCFile)
      ->required();

  sc.add_option(
      "--resolution",
      c.resolution,
      "HiC matrix resolution (ignored when file is in .cool format).")
      ->check(CLI::NonNegativeNumber)
      ->transform(AsGenomicDistance);

  sc.add_option(
      "--matrix-type",
      c.matrix_type,
      "Matrix type (ignored when file is not in .hic format).")
      ->check(CLI::IsMember{{"observed", "oe", "expected"}, CLI::ignore_case})
      ->transform(ParseHiCMatrixType)
      ->default_str("observed");

  sc.add_option(
      "--matrix-unit",
      c.matrix_unit,
      "Matrix unit (ignored when file is not in .hic format).")
      ->check(CLI::IsMember{{"BP", "FRAG"}, CLI::ignore_case})
      ->transform(ParseHiCMatrixUnit)
      ->default_str("BP");

  sc.add_option(
      "-t,--table",
      c.table,
      "Name of the table to dump.\n")
      ->check(CLI::IsMember({"chroms", "bins", "pixels", "normalizations",
                             "resolutions", "cells", "weights"}))
      ->capture_default_str();

  sc.add_option(
      "-r,--range",
      c.range1,
      "Coordinates of the genomic regions to be dumped following UCSC style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--range2",
      c.range2,
      "Coordinates of the genomic regions to be dumped following UCSC style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--query-file",
      c.query_file,
      "Path to a BEDPE file with the list of coordinates to be fetched (pass - to read queries from stdin).")
      ->check(CLI::ExistingFile | CLI::IsMember({"-"}))
      ->capture_default_str();

  sc.add_flag(
      "--cis-only",
      c.cis_only,
      "Dump intra-chromosomal interactions only.");

  sc.add_flag(
      "--trans-only",
      c.trans_only,
      "Dump inter-chromosomal interactions only.");

  sc.add_option(
      "-b,--balance",
      c.normalization,
      "Balance interactions using the given method.")
      ->capture_default_str();

  sc.add_flag(
      "--sorted,!--unsorted",
      c.sorted,
      "Return interactions in ascending order.")
      ->capture_default_str();

  sc.add_flag(
      "--join,!--no-join",
      c.join,
      "Output pixels in BG2 format.")
      ->capture_default_str();
  // clang-format on

  sc.get_option("--range2")->needs(sc.get_option("--range"));

  sc.get_option("--query-file")->excludes(sc.get_option("--range"));
  sc.get_option("--query-file")->excludes(sc.get_option("--range2"));
  sc.get_option("--query-file")->excludes(sc.get_option("--cis-only"));
  sc.get_option("--query-file")->excludes(sc.get_option("--trans-only"));

  sc.get_option("--cis-only")->excludes(sc.get_option("--trans-only"));
  sc.get_option("--cis-only")->excludes(sc.get_option("--range"));
  sc.get_option("--cis-only")->excludes(sc.get_option("--range2"));

  sc.get_option("--trans-only")->excludes(sc.get_option("--range"));
  sc.get_option("--trans-only")->excludes(sc.get_option("--range2"));

  _config = std::monostate{};
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void Cli::validate_dump_subcommand() const {
  assert(_cli.get_subcommand("dump")->parsed());

  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(_config);

  const auto& subcmd = *_cli.get_subcommand("dump");

  const auto is_hic = hic::utils::is_hic_file(c.uri);
  const auto is_cooler = cooler::utils::is_cooler(c.uri);
  const auto is_mcooler = cooler::utils::is_multires_file(c.uri);
  const auto is_scool = cooler::utils::is_scool_file(c.uri);

  if ((is_hic || is_mcooler) && !c.resolution.has_value() &&
      (c.table == "pixels" || c.table == "bins" || c.table == "weights")) {
    const auto resolutions =
        is_hic ? hic::utils::list_resolutions(c.uri) : cooler::utils::list_resolutions(c.uri);
    if (resolutions.size() != 1) {
      errors.emplace_back("--resolution is mandatory when file is in .hic or .mcool format.");
    }
  }

  const auto resolution_parsed = !subcmd.get_option("--resolution")->empty();

  if ((is_cooler || is_scool) && resolution_parsed) {
    _warnings.emplace_back("--resolution is ignored when file is in .[s]cool format.");
  }

  const auto range_parsed = !subcmd.get_option("--range")->empty();
  if (range_parsed && c.table != "chroms" && c.table != "bins" && c.table != "pixels" &&
      c.table != "weights") {
    _warnings.emplace_back(
        "--range and --range2 are ignored when --table is not bins, chroms, pixels, or weights");
  }

  const auto query_file_parsed = !subcmd.get_option("--query-file")->empty();
  if (query_file_parsed && c.table != "bins" && c.table != "pixels") {
    _warnings.emplace_back("--query-file is ignored when --table is not bins or pixels");
  }

  const auto matrix_type_parsed = !subcmd.get_option("--matrix-type")->empty();
  const auto matrix_unit_parsed = !subcmd.get_option("--matrix-unit")->empty();

  if (!is_hic && (matrix_type_parsed || matrix_unit_parsed)) {
    _warnings.emplace_back(
        "--matrix-type and --matrix-unit are ignored when input file is not in .hic format.");
  }

  if (is_hic && c.matrix_unit == hic::MatrixUnit::FRAG) {
    errors.emplace_back("--matrix-type=FRAG is not yet supported.");
  }

  if ((c.cis_only || c.trans_only) && c.table != "pixels") {
    errors.emplace_back("--cis-only and --trans-only require --table=pixels.");
  }

  const auto join_parsed = !subcmd.get_option("--join")->empty();
  if (join_parsed && c.table != "pixels") {
    errors.emplace_back("--join requires --table=pixels.");
  }

  const auto balance_parsed = !subcmd.get_option("--balance")->empty();
  if (balance_parsed && c.table != "pixels") {
    errors.emplace_back("--balance requires --table=pixels.");
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}\n"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_dump_subcommand() {
  auto& c = std::get<DumpConfig>(_config);

  c.format = infer_input_format(c.uri);
  if (c.format == "scool" && (c.table == "chroms" || c.table == "bins")) {
    const cooler::SingleCellFile sclr{c.uri};
    if (sclr.cells().empty()) {
      throw std::runtime_error("file does not contain any cell");
    }
    c.uri = fmt::format(FMT_STRING("{}::/cells/{}"), c.uri, *sclr.cells().begin());
    c.format = "cool";
  }

  if (c.table != "bins" && c.table != "pixels") {
    c.query_file = "";
  }

  if (_cli.get_subcommand("dump")->get_option("--range2")->empty()) {
    c.range2 = c.range1;
  }

  // in spdlog, high numbers correspond to low log levels
  assert(c.verbosity > 0 && c.verbosity < 5);  // NOLINTNEXTLINE(*-narrowing-conversions)
  c.verbosity = parse_hictk_verbosity_from_env().value_or(
      static_cast<std::int16_t>(spdlog::level::critical) - c.verbosity);
}

}  // namespace hictk::tools
