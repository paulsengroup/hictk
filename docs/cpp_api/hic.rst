..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

.. cpp:namespace:: hictk

Hi-C API
##########

API to operate on .hic files. Compared to the generic API, this API provides:

* more control over how files are opened
* access to .hic-specific metadata
* control over the interaction block cache

Common
------

.. cpp:namespace:: hictk::hic

.. cpp:enum-class:: MatrixType

  .. cpp:enumerator:: observed

  .. cpp:enumerator:: oe

  .. cpp:enumerator:: expected

.. cpp:enum-class:: MatrixUnit

  .. cpp:enumerator:: BP

  .. cpp:enumerator:: FRAG

.. cpp:enum-class:: QUERY_TYPE

  .. cpp:enumerator:: BED

  .. cpp:enumerator:: UCSC


File handle
-----------

.. cpp:namespace:: hictk::hic

.. cpp:class:: File

  **Constructors**

  .. cpp:function:: explicit File(std::string url_, std::uint32_t resolution_, MatrixType type_ = MatrixType::observed, MatrixUnit unit_ = MatrixUnit::BP, std::uint64_t block_cache_capacity = 0);

  **Open/close methods**

  .. cpp:function:: File &open(std::string url_, std::uint32_t resolution_, MatrixType type_ = MatrixType::observed, MatrixUnit unit_ = MatrixUnit::BP, std::uint64_t block_cache_capacity = 0);
  .. cpp:function:: File &open(std::uint32_t resolution_, MatrixType type_ = MatrixType::observed, MatrixUnit unit_ = MatrixUnit::BP, std::uint64_t block_cache_capacity = 0);

  **Accessors**

  .. cpp:function:: [[nodiscard]] bool has_resolution(std::uint32_t resolution) const;

  .. cpp:function:: [[nodiscard]] const std::string &path() const noexcept;
  .. cpp:function:: [[nodiscard]] const std::string &name() const noexcept;

  .. cpp:function:: [[nodiscard]] std::int32_t version() const noexcept;

  .. cpp:function:: [[nodiscard]] const Reference &chromosomes() const noexcept;
  .. cpp:function:: [[nodiscard]] const BinTable &bins() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;

  .. cpp:function:: [[nodiscard]] std::uint32_t resolution() const noexcept;
  .. cpp:function:: [[nodiscard]] std::uint64_t nbins() const;
  .. cpp:function:: [[nodiscard]] std::uint64_t nchroms(bool include_ALL = false) const;
  .. cpp:function:: [[nodiscard]] const std::string &assembly() const noexcept;
  .. cpp:function:: [[nodiscard]] const std::vector<std::uint32_t> &avail_resolutions() const noexcept;
  .. cpp:function:: [[nodiscard]] bool has_normalization(std::string_view normalization) const;
  .. cpp:function:: [[nodiscard]] std::vector<balancing::Method> avail_normalizations() const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(const balancing::Method &norm, const Chromosome &chrom) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(std::string_view norm, const Chromosome &chrom) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(const balancing::Method &norm) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(std::string_view norm) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(const balancing::Method &norm, const Chromosome &chrom) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(std::string_view norm, const Chromosome &chrom) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(const balancing::Method &norm) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(std::string_view norm) const;

  **Fetch methods (1D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelectorAll fetch(const balancing::Method &norm = balancing::Method::NONE()) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range, const balancing::Method &norm = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start, std::uint32_t end, const balancing::Method &norm = balancing::Method::NONE()) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::uint64_t first_bin, std::uint64_t last_bin, const balancing::Method &norm = balancing::Method::NONE()) const;

  **Fetch methods (2D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2, const balancing::Method &norm = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1, std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2, const balancing::Method &norm = balancing::Method::NONE()) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::uint64_t first_bin1, std::uint64_t last_bin1, std::uint64_t first_bin2, std::uint64_t last_bin2, const balancing::Method &norm = balancing::Method::NONE()) const;

  **Caching**

  .. cpp:function:: [[nodiscard]] std::size_t num_cached_footers() const noexcept;
  .. cpp:function:: void purge_footer_cache();

  .. cpp:function:: [[nodiscard]] double block_cache_hit_rate() const noexcept;
  .. cpp:function:: void reset_cache_stats() const noexcept;
  .. cpp:function:: void clear_cache() noexcept;
  .. cpp:function:: void optimize_cache_size(std::size_t upper_bound = (std::numeric_limits<std::size_t>::max)());
  .. cpp:function:: void optimize_cache_size_for_iteration(std::size_t upper_bound = (std::numeric_limits<std::size_t>::max)());
  .. cpp:function:: void optimize_cache_size_for_random_access(std::size_t upper_bound = (std::numeric_limits<std::size_t>::max)());
  .. cpp:function:: [[nodiscard]] std::size_t cache_capacity() const noexcept;

Pixel selector
--------------

.. cpp:namespace:: hictk::hic

.. cpp:class:: PixelSelector

  **Operators**

  .. cpp:function:: [[nodiscard]] bool operator==(const PixelSelector &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const PixelSelector &other) const noexcept;

  **Iteration**

  .. cpp:function:: template <typename N> [[nodiscard]] auto begin(bool sorted = true) const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto end() const -> iterator<N>;

  .. cpp:function:: template <typename N> [[nodiscard]] auto cbegin(bool sorted = true) const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto cend() const -> iterator<N>;

  **Fetch at once**

  .. cpp:function:: template <typename N> [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  **Accessors**

  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  .. cpp:function:: [[nodiscard]] MatrixType matrix_type() const noexcept;
  .. cpp:function:: [[nodiscard]] const balancing::Method& normalization() const noexcept;
  .. cpp:function:: [[nodiscard]] MatrixUnit unit() const noexcept;
  .. cpp:function:: [[nodiscard]] std::uint32_t resolution() const noexcept;

  .. cpp:function:: [[nodiscard]] const Chromosome &chrom1() const noexcept;
  .. cpp:function:: [[nodiscard]] const Chromosome &chrom2() const noexcept;

  .. cpp:function:: [[nodiscard]] const balancing::Weights &weights1() const noexcept;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &weights2() const noexcept;

  .. cpp:function:: [[nodiscard]] const BinTable &bins() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
  .. cpp:function:: [[nodiscard]] const internal::HiCFooterMetadata &metadata() const noexcept;

  .. cpp:function:: [[nodiscard]] bool is_inter() const noexcept;
  .. cpp:function:: [[nodiscard]] bool is_intra() const noexcept;
  .. cpp:function:: [[nodiscard]] bool empty() const noexcept;

  **Caching**

  .. cpp:function:: [[nodiscard]] std::size_t estimate_optimal_cache_size(std::size_t num_samples = 500) const;
  .. cpp:function:: void clear_cache() const;
