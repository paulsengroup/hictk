// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

[[nodiscard]] static phmap::btree_map<std::string, std::string>
generate_mappings_add_chr_prefix_prefix(std::string_view uri) {
  phmap::btree_map<std::string, std::string> mappings{};
  const auto chroms = cooler::File{uri}.chromosomes();
  for (const auto& chrom : chroms) {
    mappings.emplace(std::string{chrom.name()}, "chr" + std::string{chrom.name()});
  }
  return mappings;
}

[[nodiscard]] static phmap::btree_map<std::string, std::string>
generate_mappings_remove_chr_prefix_prefix(std::string_view uri) {
  phmap::btree_map<std::string, std::string> mappings{};
  const auto chroms = cooler::File{uri}.chromosomes();
  for (const auto& chrom : chroms) {
    const auto match = chrom.name().find("chr") == 0;
    if (match) {
      mappings.emplace(std::string{chrom.name()}, std::string{chrom.name().substr(3)});
    }
  }
  return mappings;
}

[[nodiscard]] static phmap::btree_map<std::string, std::string> read_mappings_from_file(
    const std::filesystem::path& path) {
  if (path.empty()) {
    return {};
  }

  std::ifstream ifs;
  ifs.exceptions(std::ios::badbit);
  ifs.open(path);

  phmap::btree_map<std::string, std::string> mappings{};
  std::string buff{};

  for (std::size_t i = 0; std::getline(ifs, buff); ++i) {
    if (buff.empty()) {
      continue;
    }
    if (buff.back() == '\r') {
      buff = buff.substr(0, buff.size() - 1);
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

[[nodiscard]] static phmap::btree_map<std::string, std::string> generate_name_mappings(
    std::string_view uri, const std::filesystem::path& name_mappings_path, bool add_chr_prefix,
    [[maybe_unused]] bool remove_chr_prefix) {
  phmap::btree_map<std::string, std::string> mappings{};
  if (!name_mappings_path.empty()) {
    mappings = read_mappings_from_file(name_mappings_path);
  }
  if (add_chr_prefix) {
    mappings = generate_mappings_add_chr_prefix_prefix(uri);
  }

  if (remove_chr_prefix) {
    mappings = generate_mappings_remove_chr_prefix_prefix(uri);
  }

  [[maybe_unused]] std::string mappings_str{};
  std::for_each(mappings.begin(), mappings.end(), [&](const auto& m) {
    mappings_str += fmt::format(FMT_STRING("\n - {} -> {}"), m.first, m.second);
  });

  if (mappings.empty()) {
    SPDLOG_WARN("Chromosome name map is empty: no chromosomes will be renamed!");
  } else {
    SPDLOG_INFO(FMT_STRING("Renaming chromosomes as follows:{}"), mappings_str);
  }

  return mappings;
}

static void remove_hardlinks_scool(HighFive::File& h5f,
                                   const phmap::btree_set<std::string>& cells) {
  for (const auto& cell : cells) {
    h5f.unlink(fmt::format(FMT_STRING("/cells/{}/chroms"), cell));
  }
}

static void create_hardlinks_scool(HighFive::File& h5f,
                                   const phmap::btree_set<std::string>& cells) {
  const auto chrom_grp = h5f.getGroup("/chroms");
  for (const auto& cell : cells) {
    h5f.createHardLink(fmt::format(FMT_STRING("/cells/{}/chroms"), cell), chrom_grp);
  }
}

[[nodiscard]] static int rename_chromosomes_cooler(const RenameChromosomesConfig& c) {
  const auto mappings =
      generate_name_mappings(c.uri, c.path_to_name_mappings, c.add_chr_prefix, c.remove_chr_prefix);
  cooler::utils::rename_chromosomes(c.uri, mappings);
  return 0;
}

[[nodiscard]] static int rename_chromosomes_multires_cooler(const RenameChromosomesConfig& c) {
  const auto resolutions = cooler::MultiResFile(c.uri).resolutions();
  const auto mappings = generate_name_mappings(
      fmt::format(FMT_STRING("{}::/resolutions/{}"), c.uri, resolutions.front()),
      c.path_to_name_mappings, c.add_chr_prefix, c.remove_chr_prefix);
  for (const auto& res : resolutions) {
    cooler::utils::rename_chromosomes(fmt::format(FMT_STRING("{}::/resolutions/{}"), c.uri, res),
                                      mappings);
  }
  return 0;
}

[[nodiscard]] static int rename_chromosomes_single_cell_cooler(const RenameChromosomesConfig& c) {
  assert(cooler::utils::is_scool_file(c.uri));
  const auto cells = cooler::SingleCellFile(c.uri).cells();
  const auto uri = fmt::format(FMT_STRING("{}::/cells/{}"), c.uri, *cells.begin());

  const auto mappings =
      generate_name_mappings(uri, c.path_to_name_mappings, c.add_chr_prefix, c.remove_chr_prefix);

  // NOLINTNEXTLINE(misc-const-correctness)
  HighFive::File h5f(c.uri, HighFive::File::ReadWrite);

  remove_hardlinks_scool(h5f, cells);

  const cooler::RootGroup root_grp{h5f.getGroup("/")};
  cooler::Dataset dset{root_grp, "/chroms/name"};
  cooler::utils::rename_chromosomes(dset, mappings);

  create_hardlinks_scool(h5f, cells);
  assert(cooler::utils::is_scool_file(c.uri));

  return 0;
}

int rename_chromosomes_subcmd(const RenameChromosomesConfig& c) {
  if (cooler::utils::is_cooler(c.uri)) {
    return rename_chromosomes_cooler(c);
  }

  if (cooler::utils::is_multires_file(c.uri)) {
    return rename_chromosomes_multires_cooler(c);
  }
  return rename_chromosomes_single_cell_cooler(c);
}

}  // namespace hictk::tools
