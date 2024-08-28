// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

#include "hictk/chromosome.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/tools/common.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void update_status_table(const cooler::utils::ValidationStatusCooler& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("is_valid_cooler", status.is_cooler);
}

static void update_status_table(const cooler::utils::ValidationStatusMultiresCooler& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("is_valid_mcool", status.is_multires_file);

  assert(status.valid_resolutions.empty());
  assert(status.invalid_resolutions.empty());
}

static void update_status_table(const cooler::utils::ValidationStatusScool& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("unexpected_number_of_cells", status.unexpected_number_of_cells);
  buff.insert("is_valid_scool", status.is_scool_file);

  assert(status.valid_cells.empty());
  assert(status.invalid_cells.empty());
}

static std::optional<hic::File> open_hic_noexcept(const std::string& path,
                                                  std::uint32_t resolution) noexcept {
  try {
    return hic::File(path, resolution);
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("[{}] failed to open file: {}"), resolution, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("[{}] failed to open file: unknown error"), resolution);
    return {};
  }
}

[[nodiscard]] static std::optional<cooler::MultiResFile> open_mcool_noexcept(
    std::string_view uri) noexcept {
  try {
    return cooler::MultiResFile{std::filesystem::path(uri)};
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": {}"), uri, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": unknown error"), uri);
    return {};
  }
}

[[nodiscard]] static std::optional<cooler::SingleCellFile> open_scool_noexcept(
    std::string_view uri) noexcept {
  try {
    return cooler::SingleCellFile{std::filesystem::path(uri)};
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": {}"), uri, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": unknown error"), uri);
    return {};
  }
}

static void validate_hic(const hic::File& hf, const Chromosome& chrom1, const Chromosome& chrom2) {
  if (chrom1.is_all() || chrom2.is_all()) {
    return;
  }

  try {
    std::ignore = hf.fetch(chrom1.name(), chrom2.name());
  } catch (const std::exception& e) {
    const std::string_view msg{e.what()};
    if (msg.find("Unable to find block map") != std::string_view::npos) {
      return;
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Validation failed for {}:{} map at {} resolution: {}"),
                    chrom1.name(), chrom2.name(), hf.resolution(), e.what()));
  }
}

[[nodiscard]] static std::pair<int, toml::table> validate_hic(const std::string& path,
                                                              bool exhaustive) {
  toml::table status;
  int return_code = 0;
  for (const auto& res : hic::utils::list_resolutions(path)) {
    const auto hf = open_hic_noexcept(path, res);
    if (!hf) {
      status.insert(fmt::to_string(res), "unable to open resolution");
      return_code = 1;
      assert(false);  // This should never happen
      if (!exhaustive) {
        return std::make_pair(return_code, status);
      }
      continue;
    }

    const auto& chroms = hf->chromosomes();
    for (std::uint32_t i = 0; i < chroms.size(); ++i) {
      for (std::uint32_t j = i; j < chroms.size(); ++j) {
        const auto& chrom1 = chroms.at(i);
        const auto& chrom2 = chroms.at(j);

        try {
          validate_hic(*hf, chrom1, chrom2);
        } catch ([[maybe_unused]] const std::exception& e) {
          SPDLOG_DEBUG(FMT_STRING("[{}]: validation failed for {}:{} {}"), res, chrom1.name(),
                       chrom2.name(), e.what());
          status.insert(fmt::format(FMT_STRING("{}:{}_{}"), chrom1.name(), chrom2.name(), res),
                        "unable to fetch interactions");
          return_code = 1;
          if (!exhaustive) {
            return std::make_pair(return_code, status);
          }
        }
      }
    }
  }

  return std::make_pair(return_code, status);
}

[[nodiscard]] static bool validate_bin_table_shape(const cooler::File& clr) {
  const auto& chroms = clr.dataset("bins/chrom");
  const auto& starts = clr.dataset("bins/start");
  const auto& ends = clr.dataset("bins/end");

  const auto expected_num_bins = clr.bins().size();
  if (chroms.size() != expected_num_bins || starts.size() != expected_num_bins ||
      ends.size() != expected_num_bins) {
    return false;
  }
  return true;
}

[[nodiscard]] static bool validate_bins_dtypes(const cooler::File& clr) {
  try {
    std::ignore = *clr.dataset("bins/chrom").begin<std::string>();
    std::ignore = *clr.dataset("bins/start").begin<std::int32_t>();
    std::ignore = *clr.dataset("bins/end").begin<std::int32_t>();
    return true;
  } catch (...) {
    return false;
  }
}

[[nodiscard]] static std::size_t count_invalid_bins(const cooler::File& clr) {
  const auto& chroms = clr.dataset("bins/chrom");
  const auto& starts = clr.dataset("bins/start");
  const auto& ends = clr.dataset("bins/end");

  auto first_chrom_id = chroms.begin<std::int32_t>();
  const auto last_chrom = chroms.end<std::int32_t>();
  auto first_start = starts.begin<std::int32_t>();
  const auto last_start = starts.end<std::int32_t>();
  auto first_end = ends.begin<std::int32_t>();
  const auto last_end = ends.end<std::int32_t>();

  auto first_bin = clr.bins().begin();
  const auto last_bin = clr.bins().end();

  std::size_t num_invalid_bins{};
  const auto num_bins = clr.bins().size();
  for (std::size_t i = 0; i < num_bins; ++i) {
    const auto bin = *first_bin++;

    const auto chrom_it = clr.chromosomes().find(static_cast<std::uint32_t>(*first_chrom_id++));
    if (chrom_it == clr.chromosomes().end()) {
      ++first_start;
      ++first_end;
      ++num_invalid_bins;
      continue;
    }

    const auto& chrom = *chrom_it;
    const auto start = static_cast<std::uint32_t>(*first_start++);
    const auto end = static_cast<std::uint32_t>(*first_end++);

    if (bin.chrom() != chrom || bin.start() != start || bin.end() != end) {
      ++num_invalid_bins;
    }
  }

  return num_invalid_bins;
}

static bool check_bin_table(const cooler::File& clr, toml::table& status) {
  const auto size_ok = validate_bin_table_shape(clr);
  status.insert("bin_table_shape_ok", size_ok);
  if (!size_ok) {
    return false;
  }

  const auto dtypes_ok = validate_bins_dtypes(clr);
  status.insert("bin_table_dtypes_ok", dtypes_ok);
  if (!dtypes_ok) {
    return false;
  }

  const auto num_invalid_bins = count_invalid_bins(clr);
  status.insert("bin_table_num_invalid_bins", static_cast<std::int64_t>(num_invalid_bins));
  return num_invalid_bins == 0;
}

[[nodiscard]] static std::pair<int, toml::table> validate_cooler(std::string_view path,
                                                                 bool validate_index) {
  int return_code = 0;
  toml::table status;

  update_status_table(cooler::utils::is_cooler(path), status);
  auto is_cooler = *status.get("is_valid_cooler")->as_boolean();

  std::optional<cooler::File> clr{};
  if (is_cooler) {
    try {
      clr = cooler::File(path);
    } catch (...) {
      is_cooler = false;
      status.insert_or_assign("is_cooler", is_cooler);
      return_code = 1;
    }
  }

  if (!!clr) {
    assert(clr.has_value());
    check_bin_table(*clr, status);
  } else {
    status.insert("bin_table_shape_ok", "not_checked");
    status.insert("bin_table_dtypes_ok", "not_checked");
    status.insert("bin_table_num_invalid_bins", "not_checked");
  }

  if (!!clr && validate_index) {
    std::string buff{};
    const auto index_ok = cooler::utils::index_is_valid(path, buff);
    if (index_ok) {
      status.insert("index_is_valid", true);
    } else {
      assert(!buff.empty());
      return_code = 1;
      status.insert("index_is_valid", buff);
    }
  } else {
    status.insert("index_is_valid", "not_checked");
  }

  return std::make_pair(return_code, status);
}

[[nodiscard]] static std::pair<int, toml::table> validate_mcool(std::string_view path,
                                                                bool validate_index,
                                                                bool exhaustive) {
  int return_code = 0;
  toml::table global_status;
  const auto validation_status = cooler::utils::is_multires_file(path, false);
  update_status_table(validation_status, global_status);

  if (!validation_status.is_multires_file) {
    return std::make_pair(1, global_status);
  }

  const auto mclr = open_mcool_noexcept(path);
  if (!mclr) {
    global_status.insert_or_assign("is_mcool", false);
    return std::make_pair(1, global_status);
  }

  bool early_return = false;
  std::for_each(mclr->resolutions().rbegin(), mclr->resolutions().rend(), [&](const auto res) {
    if (early_return) {
      return;
    }
    const auto [_, status] = validate_cooler(mclr->open(res).uri(), validate_index);
    const auto is_cooler = status.get("is_cooler")->as_boolean();
    global_status.insert(fmt::to_string(res), status);
    if (!is_cooler && !exhaustive) {
      return_code = 1;
      early_return = true;
    }
  });

  return std::make_pair(return_code, global_status);
}

[[nodiscard]] static std::pair<int, toml::table> validate_scool(std::string_view path,
                                                                bool validate_index,
                                                                bool exhaustive) {
  int return_code = 0;
  toml::table global_status;

  update_status_table(cooler::utils::is_scool_file(path, false), global_status);
  const auto sclr = open_scool_noexcept(path);
  if (!sclr) {
    global_status.insert_or_assign("is_valid_scool", false);
    return std::make_pair(1, global_status);
  }

  for (const auto& cell : sclr->cells()) {
    const auto [_, status] = validate_cooler(sclr->open(cell).uri(), validate_index);
    const auto is_cooler = status.get("is_cooler")->as_boolean();
    global_status.insert(cell, status);

    if (!is_cooler) {
      return_code = 1;
      if (!exhaustive) {
        std::make_pair(1, global_status);
      }
    }
  }

  return std::make_pair(return_code, global_status);
}

static void print_report(const toml::table& status, std::string_view format) {
  if (format == "json") {
    fmt::print(FMT_STRING("{}\n"), io::toml::format_to_json(status, {}));
    return;
  }

  if (format == "toml") {
    fmt::print(FMT_STRING("{}\n"), io::toml::format_to_toml(status, {}));
    return;
  }

  assert(format == "yaml");
  fmt::print(FMT_STRING("{}\n"), io::toml::format_to_yaml(status, {}));
}

[[nodiscard]] static toml::table merge_tables(const toml::table& t1, const toml::table& t2) {
  auto t = t1;
  for (const auto& [k, v] : t2) {
    t.insert(k, v);
  }
  return t;
}

int validate_subcmd(const ValidateConfig& c) {
  try {
    int return_code = 0;
    toml::table status;

    if (c.quiet) {
      // In theory nothing should write to stdout, but better to be safe than sorry
#ifdef _WIN32
      std::freopen("nul", "w", stdout);
#else
      std::freopen("/dev/null", "w", stdout);
#endif
    }

    const auto is_cooler = cooler::utils::is_cooler(c.uri);
    const auto is_hic = hic::utils::is_hic_file(c.uri);
    const auto is_mcool = cooler::utils::is_multires_file(c.uri, false);
    const auto is_scool = cooler::utils::is_scool_file(c.uri, false);

    if (c.include_file_path) {
      status.insert("uri", c.uri);
    }

    if (is_cooler) {
      status.insert("format", "cool");
    }
    if (is_hic) {
      status.insert("format", "hic");
    }
    if (is_mcool) {
      status.insert("format", "mcool");
    }
    if (is_scool) {
      status.insert("format", "scool");
    }

    if (!is_hic && !is_cooler && !is_mcool && !is_scool) {
      if (!c.quiet) {
        print_report(status, c.output_format);
        fmt::print(stderr, FMT_STRING("### FAILURE: \"{}\" is not in .hic or .[ms]cool format!\n"),
                   c.uri);
      }
      return 1;
    }

    if (is_hic) {
      auto res = validate_hic(c.uri, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else if (is_mcool) {
      auto res = validate_mcool(c.uri, c.validate_index, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else if (is_scool) {
      auto res = validate_scool(c.uri, c.validate_index, c.exhaustive);
      return_code = res.first;
      status = merge_tables(status, res.second);
    } else {
      auto res = validate_cooler(c.uri, c.validate_index);
      return_code = res.first;
      status = merge_tables(status, res.second);
    }

    if (!c.quiet) {
      print_report(status, c.output_format);
      if (is_hic) {
        fmt::print(stderr, FMT_STRING("### SUCCESS: \"{}\" is a valid .hic file."), c.uri);
      } else if (is_mcool) {
        fmt::print(stderr, FMT_STRING("### SUCCESS: \"{}\" is a valid .mcool file."), c.uri);
      } else if (is_scool) {
        fmt::print(stderr, FMT_STRING("### SUCCESS: \"{}\" is a valid .scool file."), c.uri);
      } else if (std::filesystem::exists(c.uri)) {
        fmt::print(stderr, FMT_STRING("### SUCCESS: \"{}\" is a valid .cool file."), c.uri);
      } else {
        fmt::print(stderr, FMT_STRING("### SUCCESS: \"{}\" points to valid Cooler."), c.uri);
      }
    }

    return return_code;
  } catch (...) {
    if (c.quiet) {
      return 1;
    }
    throw;
  }
}

}  // namespace hictk::tools
