// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>

#include <CLI/CLI.hpp>
#include <cassert>
#include <string>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void Cli::make_dump_subcommand() {
  auto& sc = *this->_cli.add_subcommand("dump", "Dump Cooler data to stdout.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] std::size_t i) {
                    assert(this->_config.index() == 0);
                    this->_config = DumpConfig{};
                  });

  this->_config = DumpConfig{};
  auto& c = std::get<DumpConfig>(this->_config);

  // clang-format off
  sc.add_option(
      "uri",
      c.uri,
      "Path to a .hic, .cool or .mcool file (Cooler URI syntax supported).")
      ->check(IsValidHiCFile | IsValidCoolerFile)
      ->required();

  sc.add_option(
      "--resolution",
      c.resolution,
      "HiC matrix resolution (ignored when file is not in .hic format).")
      ->check(CLI::NonNegativeNumber);

  sc.add_option(
      "--matrix-type",
      c.matrix_type,
      "Matrix type (ignored when file is not in .hic format).")
      ->transform(ParseHiCMatrixType)
      ->default_str("observed");

  sc.add_option(
      "--matrix-unit",
      c.matrix_unit,
      "Matrix unit (ignored when file is not in .hic format).")
      ->transform(ParseHiCMatrixUnit)
      ->default_str("BP");

  sc.add_option(
      "-t,--table",
      c.table,
      "Name of the table to dump.\n")
      ->check(CLI::IsMember({"chroms", "bins", "pixels"}))
      ->capture_default_str();

  sc.add_option(
      "-r,--range",
      c.range1,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--range2",
      c.range2,
      "Coordinates of the genomic regions to be dumped following UCSC-style notation (chr1:0-1000).")
      ->capture_default_str();

  sc.add_option(
      "--query-file",
      c.query_file,
      "Path to a BEDPE file with the list of coordinates to be fetched (pass - to read queries from stdin).")
      ->check(CLI::ExistingFile | CLI::IsMember({"-"}))
      ->capture_default_str();

  sc.add_option(
      "-b,--balance",
      c.normalization,
      "Balance interactions using the given method.")
      ->capture_default_str();

  sc.add_flag(
      "--join,!--no-join",
      c.join,
      "Output pixels in BG2 format.")
      ->capture_default_str();

  sc.add_option(
      "--weight-type",
      c.weight_type,
      "How balancing weights should be applied to raw interactions (ignored when file is in .hic format).")
      ->check(CLI::IsMember({"infer", "divisive", "multiplicative"}))
      ->capture_default_str();

  // clang-format on

  sc.get_option("--query-file")->excludes(sc.get_option("--range"));
  sc.get_option("--query-file")->excludes(sc.get_option("--range2"));

  this->_config = std::monostate{};
}

void Cli::validate_dump_subcommand() const {
  assert(this->_cli.get_subcommand("dump")->parsed());

  [[maybe_unused]] std::vector<std::string> warnings;
  std::vector<std::string> errors;
  const auto& c = std::get<DumpConfig>(this->_config);

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }

  const auto is_hic = hic::utils::is_hic_file(c.uri);
  const auto is_cooler = cooler::utils::is_cooler(c.uri);
  const auto is_mcooler = cooler::utils::is_multires_file(c.uri);

  if (is_hic && c.resolution == 0 && c.table != "chroms") {
    errors.emplace_back("--resolution is mandatory when file is in .hic format.");
  }

  const auto resolution_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--resolution")->empty();

  if ((is_cooler || is_mcooler) && resolution_parsed) {
    warnings.emplace_back("--resolution is ignored when file is in .cool or .mcool format.");
  }

  const auto weight_type_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--weight-type")->empty();

  if (is_hic && weight_type_parsed) {
    warnings.emplace_back("--weight-type is ignored when file is in .hic format.");
  }

  const auto matrix_type_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--matrix-type")->empty();
  const auto matrix_unit_parsed =
      !this->_cli.get_subcommand("dump")->get_option("--matrix-unit")->empty();

  if (!is_hic && (matrix_type_parsed || matrix_unit_parsed)) {
    warnings.emplace_back(
        "--matrix-type and --matrix-unit are ignored when input file is not in .hic format.");
  }

  if (is_hic && c.matrix_unit == hic::MatrixUnit::FRAG) {
    errors.emplace_back("--matrix-type=FRAG is not yet supported.");
  }

  for (const auto& w : warnings) {
    spdlog::warn(FMT_STRING("{}"), w);
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("the following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}\n"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::transform_args_dump_subcommand() {
  auto& c = std::get<DumpConfig>(this->_config);

  c.format = infer_input_format(c.uri);
  if (c.format == "hic" && c.resolution == 0) {
    assert(c.table == "chroms");
    c.resolution = hic::utils::list_resolutions(c.uri).back();
  }

  if (this->_cli.get_subcommand("dump")->get_option("--range2")->empty()) {
    c.range2 = c.range1;
  }
}

}  // namespace hictk::tools
