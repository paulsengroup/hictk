..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

.. cpp:namespace:: hictk

Cooler API
##########

API to operate on .cool files. Compared to the generic API, this API provides:

* more control over how files are opened
* direct access to HDF5 group and datasets
* lower overhead
* support for creating .cool files
* support for opening collections of Coolers (e.g. .mcool and .scool files)

Single-resolution Cooler (.cool)
--------------------------------

.. cpp:namespace:: hictk::cooler
.. cpp:class:: File

  **Constructors**

  .. cpp:function:: File(const File &other) = delete;
  .. cpp:function:: File(File &&other) noexcept = default;

  .. cpp:function:: [[nodiscard]] explicit File(std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);
  .. cpp:function:: [[nodiscard]] explicit File(RootGroup entrypoint, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);

  **Factory functions**

  .. cpp:function:: [[nodiscard]] static File open_random_access(RootGroup entrypoint, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);
  .. cpp:function:: [[nodiscard]] static File open_read_once(RootGroup entrypoint, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);
  .. cpp:function:: template <typename PixelT = DefaultPixelT> [[nodiscard]] static File create(RootGroup entrypoint, const Reference &chroms, std::uint32_t bin_size, Attributes attributes = Attributes::init<PixelT>(0), std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, std::uint32_t compression_lvl = DEFAULT_COMPRESSION_LEVEL);
  .. cpp:function:: template <typename PixelT = DefaultPixelT> [[nodiscard]] static File create(std::string_view uri, const Reference &chroms, std::uint32_t bin_size, bool overwrite_if_exists = false, Attributes attributes = Attributes::init<PixelT>(0), std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, std::uint32_t compression_lvl = DEFAULT_COMPRESSION_LEVEL);

  **Open/close methods**

  .. cpp:function:: [[nodiscard]] static File open_random_access(std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);
  .. cpp:function:: [[nodiscard]] static File open_read_once(std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4, bool validate = true);

  .. cpp:function:: void close();

  Note that :cpp:class:`File`\s are automatically closed upon destruction.

  **Operators**

  .. cpp:function:: File &operator=(const File &other) = delete;
  .. cpp:function:: File &operator=(File &&other) noexcept = default;

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;

  Return whether the :cpp:class:`File` is in a valid state and other member functions can be safely called.

  **Accessors**

  .. cpp:function:: [[nodiscard]] std::string uri() const;
  .. cpp:function:: [[nodiscard]] std::string hdf5_path() const;
  .. cpp:function:: [[nodiscard]] std::string path() const;

  .. cpp:function:: [[nodiscard]] auto chromosomes() const noexcept -> const Reference &;
  .. cpp:function:: [[nodiscard]] auto bins() const noexcept -> const BinTable &;
  .. cpp:function:: [[nodiscard]] auto bins_ptr() const noexcept -> std::shared_ptr<const BinTable>;

  .. cpp:function:: [[nodiscard]] std::uint32_t resolution() const noexcept;
  .. cpp:function:: [[nodiscard]] std::uint64_t nbins() const;
  .. cpp:function:: [[nodiscard]] std::uint64_t nchroms() const;
  .. cpp:function:: [[nodiscard]] std::uint64_t nnz() const;

  .. cpp:function:: [[nodiscard]] auto attributes() const noexcept -> const Attributes &;
  .. cpp:function:: [[nodiscard]] auto group(std::string_view group_name) -> Group &;
  .. cpp:function:: [[nodiscard]] auto dataset(std::string_view dataset_name) -> Dataset &;
  .. cpp:function:: [[nodiscard]] auto group(std::string_view group_name) const -> const Group &;
  .. cpp:function:: [[nodiscard]] auto dataset(std::string_view dataset_name) const -> const Dataset &;

  .. cpp:function:: [[nodiscard]] const NumericVariant &pixel_variant() const noexcept;
  .. cpp:function:: template <typename T> [[nodiscard]] bool has_pixel_of_type() const noexcept;

  .. cpp:function:: [[nodiscard]] bool has_signed_pixels() const noexcept;
  .. cpp:function:: [[nodiscard]] bool has_unsigned_pixels() const noexcept;
  .. cpp:function:: [[nodiscard]] bool has_integral_pixels() const noexcept;
  .. cpp:function:: [[nodiscard]] bool has_float_pixels() const noexcept;

  **Iteration**

  .. cpp:function:: template <typename N> [[nodiscard]] typename PixelSelector::iterator<N> begin(std::string_view weight_name = "NONE") const;
  .. cpp:function:: template <typename N> [[nodiscard]] typename PixelSelector::iterator<N> end(std::string_view weight_name = "NONE") const;

  .. cpp:function:: template <typename N> [[nodiscard]] typename PixelSelector::iterator<N> cbegin(std::string_view weight_name = "NONE") const;
  .. cpp:function:: template <typename N> [[nodiscard]] typename PixelSelector::iterator<N> cend(std::string_view weight_name = "NONE") const;

  **Fetch methods (1D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(const balancing::Method &normalization = balancing::Method::NONE(), bool load_index = false) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::shared_ptr<const balancing::Weights> weights, bool load_index = false) const;

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range, std::shared_ptr<const balancing::Weights> weights, QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start, std::uint32_t end, std::shared_ptr<const balancing::Weights> weights) const;

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range, const balancing::Method &normalization = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start, std::uint32_t end, const balancing::Method &normalization = balancing::Method::NONE()) const;

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::uint64_t first_bin, std::uint64_t last_bin, std::shared_ptr<const balancing::Weights> weights = nullptr) const;

  **Fetch methods (2D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2, std::shared_ptr<const balancing::Weights> weights, QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1, std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2, std::shared_ptr<const balancing::Weights> weights) const;

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2, const balancing::Method &normalization = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1, std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2, const balancing::Method &normalization = balancing::Method::NONE()) const;

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::uint64_t first_bin1, std::uint64_t last_bin1, std::uint64_t first_bin2, std::uint64_t last_bin2, std::shared_ptr<const balancing::Weights> weights = nullptr) const;

  **Write pixels**

  .. cpp:function:: template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>> void append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate = false);

  **Normalization**

  .. cpp:function:: [[nodiscard]] bool has_normalization(std::string_view normalization) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(std::string_view normalization_, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(std::string_view normalization_, balancing::Weights::Type type, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(std::string_view normalization_, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(std::string_view normalization_, balancing::Weights::Type type, bool rescale = false) const;

  .. cpp:function:: [[nodiscard]] bool has_normalization(const balancing::Method &normalization) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(const balancing::Method &normalization_, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(const balancing::Method &normalization_, balancing::Weights::Type type, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(const balancing::Method &normalization_, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(const balancing::Method &normalization_, balancing::Weights::Type type, bool rescale = false) const;
  .. cpp:function:: [[nodiscard]] std::vector<balancing::Method> avail_normalizations() const;

  .. cpp:function:: bool purge_weights(std::string_view name = "");

  .. cpp:function:: template <typename It> static void write_weights(std::string_view uri, std::string_view name, It first_weight, It last_weight, bool overwrite_if_exists = false, bool divisive = false);
  .. cpp:function:: template <typename It> void write_weights(std::string_view name, It first_weight, It last_weight, bool overwrite_if_exists = false, bool divisive = false);

  **Others**

  .. cpp:function:: void flush();
  .. cpp:function:: void validate_bins(bool full = false) const;

Multi-resolution Cooler (.mcool)
--------------------------------

.. cpp:namespace:: hictk::cooler
.. cpp:class:: MultiResFile

  **Constructors**

  .. cpp:function:: explicit MultiResFile(const std::filesystem::path& path, unsigned int mode = HighFive::File::ReadOnly);

  **Factory functions**

  .. cpp:function:: [[nodiscard]] static MultiResFile create(const std::filesystem::path& path, const Reference& chroms, bool force_overwrite = false);
  .. cpp:function:: template <typename ResolutionIt> [[nodiscard]] static MultiResFile create(const std::filesystem::path& path, const File& base, ResolutionIt first_res, ResolutionIt last_res, bool force_overwrite = false);

  **Open/close methods**

  .. cpp:function:: [[nodiscard]] File open(std::uint32_t resolution) const;

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;

  **Accessors**

  .. cpp:function:: [[nodiscard]] std::string path() const;
  .. cpp:function:: [[nodiscard]] auto chromosomes() const noexcept -> const Reference&;
  .. cpp:function:: [[nodiscard]] constexpr const std::vector<std::uint32_t>& resolutions() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr const MultiResAttributes& attributes() const noexcept;
  .. cpp::function:: [[nodiscard]] const std::vector<balancing::Method>& avail_normalizations(std::string_view policy = "union") const;

  **Modifiers**

  .. cpp:function:: File copy_resolution(const cooler::File& clr);
  .. cpp:function:: template <typename N = DefaultPixelT> File create_resolution(std::uint32_t resolution, Attributes attributes = Attributes::init<N>(0));
  .. cpp:function:: RootGroup init_resolution(std::uint32_t resolution);

  **Others**

  .. cpp:function:: [[nodiscard]] static std::uint32_t compute_base_resolution(const std::vector<std::uint32_t>& resolutions, std::uint32_t target_res);
  .. cpp:function:: template <typename N = std::int32_t> static void coarsen(const File& clr1, File& clr2, std::vector<ThinPixel<N>>& buffer);

Single-cell Cooler (.scool)
---------------------------

.. cpp:namespace:: hictk::cooler
.. cpp:class:: SingleCellFile

  **Constructors**

  .. cpp:function:: explicit SingleCellFile(const std::filesystem::path& path, unsigned int mode = HighFive::File::ReadOnly);

  **Factory functions**

  .. cpp:function:: [[nodiscard]] static SingleCellFile create(const std::filesystem::path& path, const Reference& chroms, std::uint32_t bin_size, bool force_overwrite = false);

  **Open/close functions**

  .. cpp:function:: [[nodiscard]] File open(std::string_view cell) const;

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;

  **Accessors**

  .. cpp:function:: [[nodiscard]] std::string path() const;
  .. cpp:function:: [[nodiscard]] auto chromosomes() const noexcept -> const Reference&;
  .. cpp:function:: [[nodiscard]] auto bins() const noexcept -> const BinTable&;
  .. cpp:function:: [[nodiscard]] auto bins_ptr() const noexcept -> std::shared_ptr<const BinTable>;
  .. cpp:function:: [[nodiscard]] std::uint32_t resolution() const noexcept;

  .. cpp:function:: [[nodiscard]] constexpr const phmap::btree_set<std::string>& cells() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr const SingleCellAttributes& attributes() const noexcept;

  **Modifiers**

  .. cpp:function:: template <typename N> File create_cell(std::string_view cell, Attributes attrs = Attributes::init<N>(0));

  **Others**

  .. cpp:function:: template <typename N> File aggregate(std::string_view uri, bool overwrite_if_exists = false, std::size_t chunk_size = 500'000, std::size_t update_frequency = 10'000'000) const;

Pixel selector
--------------

.. cpp:class:: PixelSelector

  **Operators**

  .. cpp:function:: [[nodiscard]] bool operator==(const PixelSelector &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const PixelSelector &other) const noexcept;

  **Iteration**

  .. cpp:function:: template <typename N> [[nodiscard]] auto begin() const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto end() const -> iterator<N>;

  .. cpp:function:: template <typename N> [[nodiscard]] auto cbegin() const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto cend() const -> iterator<N>;

  **Fetch at once**

  .. cpp:function:: template <typename N> [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  **Accessors**

  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  .. cpp:function:: [[nodiscard]] const BinTable &bins() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
