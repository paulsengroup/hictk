// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

[[nodiscard]] static phmap::flat_hash_map<std::string, std::string>
generate_mappings_add_chr_prefix_prefix(std::string_view uri) {
  const auto chroms = cooler::File{uri}.chromosomes();
  phmap::flat_hash_map<std::string, std::string> mappings(chroms.size());
  for (const auto& chrom : chroms) {
    mappings.emplace(std::string{chrom.name()}, "chr" + std::string{chrom.name()});
  }
  return mappings;
}

[[nodiscard]] static phmap::flat_hash_map<std::string, std::string>
generate_mappings_remove_chr_prefix_prefix(std::string_view uri) {
  const auto chroms = cooler::File{uri}.chromosomes();
  phmap::flat_hash_map<std::string, std::string> mappings(chroms.size());
  for (const auto& chrom : chroms) {
    const auto match = chrom.name().find("chr") == 0;
    if (match) {
      mappings.emplace(std::string{chrom.name()}, std::string{chrom.name().substr(3)});
    }
  }
  return mappings;
}

[[nodiscard]] static phmap::flat_hash_map<std::string, std::string> read_mappings_from_file(
    const std::filesystem::path& path) {
  if (path.empty()) {
    return {};
  }

  std::ifstream ifs;
  ifs.exceptions(std::ios::badbit);
  ifs.open(path);

  phmap::flat_hash_map<std::string, std::string> mappings{};
  std::string buff{};

  for (std::size_t i = 0; std::getline(ifs, buff); ++i) {
    if (buff.empty()) {
      continue;
    }
    const auto sep_pos = buff.find('\t');
    if (sep_pos == std::string::npos) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found invalid record \"{}\" in file {} at line {}"), buff, path, i));
    }
    auto old_name = buff.substr(0, sep_pos);
    auto new_name = buff.substr(sep_pos + 1);

    if (old_name.empty() || new_name.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Found invalid record \"{}\" in file {} at line {}"), buff, path, i));
    }

    mappings.emplace(std::move(old_name), std::move(new_name));
  }

  return mappings;
}

[[nodiscard]] static phmap::flat_hash_map<std::string, std::string> generate_name_mappings(
    std::string_view uri, const std::filesystem::path& name_mappings_path, bool add_chr_prefix,
    bool remove_chr_prefix) {
  if (!name_mappings_path.empty()) {
    return read_mappings_from_file(name_mappings_path);
  }
  if (add_chr_prefix) {
    return generate_mappings_add_chr_prefix_prefix(uri);
  }

  assert(remove_chr_prefix);
  return generate_mappings_remove_chr_prefix_prefix(uri);
}

int rename_chromosomes_subcmd(const RenameChromosomesConfig& c) {
  if (cooler::utils::is_cooler(c.uri)) {
    const auto mappings = generate_name_mappings(c.uri, c.path_to_name_mappings, c.add_chr_prefix,
                                                 c.remove_chr_prefix);
    cooler::utils::rename_chromosomes(c.uri, mappings);
  } else if (cooler::utils::is_multires_file(c.uri)) {
    const auto resolutions = cooler::MultiResFile(c.uri).resolutions();
    const auto mappings = generate_name_mappings(
        fmt::format(FMT_STRING("{}::/resolutions/{}"), c.uri, resolutions.front()),
        c.path_to_name_mappings, c.add_chr_prefix, c.remove_chr_prefix);
    for (const auto& res : resolutions) {
      cooler::utils::rename_chromosomes(fmt::format(FMT_STRING("{}::/resolutions/{}"), c.uri, res),
                                        mappings);
    }
  } else {
    assert(cooler::utils::is_scool_file(c.uri));
    const auto cell_id = *cooler::SingleCellFile(c.uri).cells().begin();
    const auto uri = fmt::format(FMT_STRING("{}::/cells/{}"), c.uri, cell_id);
    const auto mappings =
        generate_name_mappings(uri, c.path_to_name_mappings, c.add_chr_prefix, c.remove_chr_prefix);
    cooler::utils::rename_chromosomes(uri, mappings);
  }

  return 1;
}

}  // namespace hictk::tools
