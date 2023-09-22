// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "hictk/bin_table.hpp"

namespace hictk::balancing {

class SparseMatrixView;
class SparseMatrix {
  std::vector<std::size_t> _bin1_ids{};
  std::vector<std::size_t> _bin2_ids{};
  std::vector<double> _counts{};

  std::uint32_t _chrom_id{};  // ID of the chromosome that is being procesed
  std::vector<std::size_t> _chrom_offsets{};
  std::vector<std::size_t> _bin1_offsets{};
  mutable std::vector<double> _marg{};

  static constexpr auto _gw_id = std::numeric_limits<std::uint32_t>::max();

 public:
  SparseMatrix() = default;
  explicit SparseMatrix(const BinTable& bins, std::uint32_t chrom_id = _gw_id);

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear() noexcept;
  void shrink_to_fit() noexcept;
  void finalize();

  [[nodiscard]] const std::vector<std::size_t>& bin1_ids() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>& bin2_ids() const noexcept;
  [[nodiscard]] const std::vector<double>& counts() const noexcept;
  [[nodiscard]] const std::vector<double>& margs() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>& chrom_offsets() const noexcept;

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count);

  [[nodiscard]] SparseMatrixView subset(std::uint32_t chrom_id) const;
  [[nodiscard]] SparseMatrixView view() const;

  void serialize(std::fstream& fs, int compression_lvl = 3) const;
  void deserialize(std::fstream& fs);
};

class SparseMatrixChunkedView;
class SparseMatrixChunked {
  mutable SparseMatrix _matrix{};
  mutable std::string _buff{};
  std::filesystem::path _path{};
  mutable std::fstream _fs{};

  std::vector<std::streamoff> _index{};

  // chrom_id, <first_offset_idx, last_offset_idx>
  phmap::flat_hash_map<std::uint32_t, std::pair<std::size_t, std::size_t>> _chrom_index{};
  std::uint32_t _chrom_id{};  // id of the chromosome that is currently being processed;
  std::size_t _size{};

  mutable std::vector<double> _marg{};
  std::vector<std::size_t> _chrom_offsets{};
  std::vector<std::size_t> _bin1_offsets{};

  std::size_t _chunk_size{};
  int _compression_lvl{};

 public:
  SparseMatrixChunked() = default;
  SparseMatrixChunked(const BinTable& bins, std::filesystem::path tmp_file, std::size_t chunk_size,
                      int compression_lvl = 3);

  SparseMatrixChunked(const SparseMatrixChunked& other) = delete;
  SparseMatrixChunked(SparseMatrixChunked&& other) noexcept = default;

  ~SparseMatrixChunked() noexcept;

  SparseMatrixChunked& operator=(const SparseMatrixChunked& other) = delete;
  SparseMatrixChunked& operator=(SparseMatrixChunked&& other) noexcept = default;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] const std::vector<double>& margs() const noexcept;
  [[nodiscard]] const std::vector<std::size_t>& chrom_offsets() const noexcept;

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count);
  void finalize();
  void finalize_chromosome(std::uint32_t chrom_id);

  [[nodiscard]] SparseMatrixChunkedView subset(std::uint32_t chrom_id) const;
  [[nodiscard]] SparseMatrixChunkedView view() const;

  void read_chunk(std::size_t chunk_id, SparseMatrix& buffer);

 private:
  void write_chunk();
};

class SparseMatrixView {
  mutable std::vector<double> _marg{};
  std::size_t _bin1_offset{};

 public:
  nonstd::span<const std::size_t> bin1_ids{};  // NOLINT
  nonstd::span<const std::size_t> bin2_ids{};  // NOLINT
  nonstd::span<const double> counts{};         // NOLINT

  SparseMatrixView() = default;
  SparseMatrixView(nonstd::span<const std::size_t> bin1_ids_,
                   nonstd::span<const std::size_t> bin2_ids_, nonstd::span<const double> counts_,
                   std::size_t bin1_offset, std::size_t num_bins);

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] const std::vector<double>& margs() const noexcept;

  const std::vector<double>& marginalize() const;
  const std::vector<double>& marginalize_nnz() const;
  const std::vector<double>& times_outer_product_marg(nonstd::span<const double> biases,
                                                      nonstd::span<const double> weights) const;
};

class SparseMatrixChunkedView {
  mutable SparseMatrix _matrix{};
  mutable std::string _buff{};
  mutable std::fstream _fs{};

  std::vector<std::streamoff> _index{};

  mutable std::vector<double> _marg{};
  std::size_t _bin1_offset{};

 public:
  SparseMatrixChunkedView() = default;
  SparseMatrixChunkedView(const std::filesystem::path& path,
                          nonstd::span<const std::streamoff> index, std::size_t bin1_offset,
                          std::size_t num_bins);

  [[nodiscard]] bool empty() const noexcept;

  [[nodiscard]] const std::vector<double>& margs() const noexcept;

  const std::vector<double>& marginalize() const;
  const std::vector<double>& marginalize_nnz() const;
  const std::vector<double>& times_outer_product_marg(nonstd::span<const double> biases,
                                                      nonstd::span<const double> weights) const;
};

}  // namespace hictk::balancing

#include "./impl/sparse_matrix_impl.hpp"
