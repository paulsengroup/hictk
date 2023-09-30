..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

.. cpp:namespace:: hictk

Shared Types
############

Types documented in this page are used throughout hictk code-base to model various concepts such as genomic intervals, reference genomes, bins and pixels.

Chromosome
----------

.. cpp:namespace:: hictk

.. cpp:class:: Chromosome

  This class models chromosomes as triplets consisting of:

  * A numeric identifier
  * The chromosome name
  * The chromosome size

  :cpp:class:`Chromosome`\s are compared by ID.

  **Constructors**

  .. cpp:function:: Chromosome() = default;
  .. cpp:function:: Chromosome(std::uint32_t id_, std::string name_, std::uint32_t size_) noexcept;

  **Operators**

  .. cpp:function:: [[nodiscard]] constexpr explicit operator bool() const noexcept;

  **Accessors**

  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t id() const noexcept;
  .. cpp:function:: [[nodiscard]] std::string_view name() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t size() const noexcept;
  .. cpp:function:: [[nodiscard]] bool is_all() const noexcept;

  **Comparison operators**

  .. cpp:function:: [[nodiscard]] constexpr bool operator<(const Chromosome& other) const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr bool operator>(const Chromosome& other) const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr bool operator<=(const Chromosome& other) const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr bool operator>=(const Chromosome& other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator==(const Chromosome& other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const Chromosome& other) const noexcept;

  .. cpp:function:: friend bool operator==(const Chromosome& a, std::string_view b_name) noexcept;
  .. cpp:function:: friend bool operator!=(const Chromosome& a, std::string_view b_name) noexcept;

  .. cpp:function:: friend bool operator==(std::string_view a_name, const Chromosome& b) noexcept;
  .. cpp:function:: friend bool operator!=(std::string_view a_name, const Chromosome& b) noexcept;

  .. cpp:function:: friend constexpr bool operator<(const Chromosome& a, std::uint32_t b_id) noexcept;
  .. cpp:function:: friend constexpr bool operator>(const Chromosome& a, std::uint32_t b_id) noexcept;
  .. cpp:function:: friend constexpr bool operator<=(const Chromosome& a, std::uint32_t b_id) noexcept;
  .. cpp:function:: friend constexpr bool operator>=(const Chromosome& a, std::uint32_t b_id) noexcept;
  .. cpp:function:: friend constexpr bool operator==(const Chromosome& a, std::uint32_t b_id) noexcept;
  .. cpp:function:: friend constexpr bool operator!=(const Chromosome& a, std::uint32_t b_id) noexcept;

  .. cpp:function:: friend constexpr bool operator<(std::uint32_t a_id, const Chromosome& b) noexcept;
  .. cpp:function:: friend constexpr bool operator>(std::uint32_t a_id, const Chromosome& b) noexcept;
  .. cpp:function:: friend constexpr bool operator<=(std::uint32_t a_id, const Chromosome& b) noexcept;
  .. cpp:function:: friend constexpr bool operator>=(std::uint32_t a_id, const Chromosome& b) noexcept;
  .. cpp:function:: friend constexpr bool operator==(std::uint32_t a_id, const Chromosome& b) noexcept;
  .. cpp:function:: friend constexpr bool operator!=(std::uint32_t a_id, const Chromosome& b) noexcept;


Genomic intervals
-----------------

.. cpp:namespace:: hictk

.. cpp:class:: GenomicInterval

  Class to represent 1D genomic intervals.

  This class has two main purposes:

  * Storing information regarding genomic intervals
  * Simplifying comparison of genomic intervals (e.g. is interval A upstream of interval B)

  .. cpp:enum-class:: QUERY_TYPE

    .. cpp:enumerator:: BED
    .. cpp:enumerator:: UCSC

  **Constructors**

  .. cpp:function:: constexpr GenomicInterval() = default;
  .. cpp:function:: explicit GenomicInterval(const Chromosome &chrom_) noexcept;
  .. cpp:function:: GenomicInterval(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end) noexcept;

  **Factory methods**

  .. cpp:function:: [[nodiscard]] static GenomicInterval parse(const Reference &chroms, std::string query, Type type = Type::UCSC);
  .. cpp:function:: [[nodiscard]] static GenomicInterval parse_ucsc(const Reference &chroms, std::string query);
  .. cpp:function:: [[nodiscard]] static GenomicInterval parse_bed(const Reference &chroms, std::string_view query, char sep = '\t');

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator==(const GenomicInterval &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const GenomicInterval &other) const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator<(const GenomicInterval &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<=(const GenomicInterval &other) const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator>(const GenomicInterval &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>=(const GenomicInterval &other) const noexcept;

  **Accessors**

  .. cpp:function:: [[nodiscard]] const Chromosome &chrom() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t start() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t end() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t size() const noexcept;


Genomic bins
------------

.. cpp:namespace:: hictk

.. cpp:class:: Bin

  Class modeling genomic bins.

  The class is implemented as a thin wrapper around :cpp:class:`GenomicInterval`\s. The main difference between :cpp:class:`Bin` and :cpp:class:`GenomicInterval` objects is that in addition to genomic coordinates, the :cpp:class:`Bin` object also store two identifiers:

  * A unique identifier that can be used to refer :cpp:class:`Bin`\s in a :cpp:class:`Reference`.
  * A relative identifier that can be uset to refer to :cpp:class:`Bin`\s in a :cpp:class:`Chromosome`.

  .. cpp:function:: constexpr Bin() = default;
  .. cpp:function:: Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  .. cpp:function:: Bin(std::uint64_t id_, std::uint32_t rel_id_, const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end_) noexcept;
  .. cpp:function:: explicit Bin(GenomicInterval interval) noexcept;
  .. cpp:function:: Bin(std::uint64_t id_, std::uint32_t rel_id_, GenomicInterval interval) noexcept;

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator==(const Bin &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const Bin &other) const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator<(const Bin &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<=(const Bin &other) const noexcept;

  .. cpp:function:: [[nodiscard]] bool operator>(const Bin &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>=(const Bin &other) const noexcept;

  .. cpp:function:: [[nodiscard]] constexpr std::uint64_t id() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t rel_id() const noexcept;
  .. cpp:function:: [[nodiscard]] const GenomicInterval &interval() const noexcept;
  .. cpp:function:: [[nodiscard]] const Chromosome &chrom() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t start() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t end() const noexcept;

  .. cpp:function:: [[nodiscard]] constexpr bool has_null_id() const noexcept;


Reference genome
----------------

.. cpp:namespace:: hictk

.. cpp:class:: Reference

  This class models the reference genome used as coordinate system in Hi-C matrices.

  :cpp:class:`Reference` objects consist of collections of :cpp:class:`Chromosome`\s with unique IDs.

  :cpp:class:`Chromosome`\s can be queried by ID or by name.

  As a general rule, queries by :cpp:class:`Chromosome` ID are more efficient than queries by name.

  **Constructors**

  .. cpp:function:: Reference() = default;

  .. cpp:function:: template <typename ChromosomeNameIt, typename ChromosomeSizeIt> Reference(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name, ChromosomeSizeIt first_chrom_size);
  .. cpp:function:: template <typename ChromosomeIt> Reference(ChromosomeIt first_chrom, ChromosomeIt last_chrom);
  .. cpp:function:: Reference(std::initializer_list<Chromosome> chromosomes);

  **Factory methods**

  .. cpp:function:: [[nodiscard]] static Reference from_chrom_sizes(const std::filesystem::path& path_to_chrom_sizes);

  **Operators**

  .. cpp:function:: [[nodiscard]] bool operator==(const Reference& other) const;
  .. cpp:function:: [[nodiscard]] bool operator!=(const Reference& other) const;

  **Iteration**

  .. cpp:function:: [[nodiscard]] auto begin() const -> const_iterator;
  .. cpp:function:: [[nodiscard]] auto end() const -> const_iterator;
  .. cpp:function:: [[nodiscard]] auto cbegin() const -> const_iterator;
  .. cpp:function:: [[nodiscard]] auto cend() const -> const_iterator;

  .. cpp:function:: [[nodiscard]] auto rbegin() const -> const_reverse_iterator;
  .. cpp:function:: [[nodiscard]] auto rend() const -> const_reverse_iterator;
  .. cpp:function:: [[nodiscard]] auto rcbegin() const -> const_reverse_iterator;
  .. cpp:function:: [[nodiscard]] auto rcend() const -> const_reverse_iterator;

  **Accessors**

  .. cpp:function:: [[nodiscard]] bool empty() const noexcept;
  .. cpp:function:: [[nodiscard]] std::size_t size() const noexcept;

  **Lookup**

  .. cpp:function:: [[nodiscard]] auto find(std::uint32_t id) const -> const_iterator;
  .. cpp:function:: [[nodiscard]] auto find(std::string_view chrom_name) const -> const_iterator;
  .. cpp:function:: [[nodiscard]] auto find(const Chromosome& chrom) const -> const_iterator;

  .. cpp:function:: [[nodiscard]] const Chromosome& at(std::uint32_t id) const;
  .. cpp:function:: [[nodiscard]] const Chromosome& at(std::string_view chrom_name) const;

  .. cpp:function:: [[nodiscard]] const Chromosome& operator[](std::uint32_t id) const noexcept;
  .. cpp:function:: [[nodiscard]] const Chromosome& operator[](std::string_view chrom_name) const noexcept;

  .. cpp:function:: [[nodiscard]] bool contains(std::uint32_t id) const;
  .. cpp:function:: [[nodiscard]] bool contains(const Chromosome& chrom) const;
  .. cpp:function:: [[nodiscard]] bool contains(std::string_view chrom_name) const;

  .. cpp:function:: [[nodiscard]] std::uint32_t get_id(std::string_view chrom_name) const;

  .. cpp:function:: [[nodiscard]] const Chromosome& longest_chromosome() const;
  .. cpp:function:: [[nodiscard]] const Chromosome& chromosome_with_longest_name() const;


Bin Table
---------

.. cpp:namespace:: hictk

.. cpp:class:: BinTable

  This class models the bin table used as coordinate system in Hi-C matrices.

  The class API gives the illusion of operating over a collection of :cpp:class:`Bin`\s.
  In reality :cpp:class:`BinTable`\s do not store any :cpp:class:`Bin`\s. All queries are satisfied through simple arithmetic operations on the prefix sum of :cpp:class:`Chromosome` sizes and :cpp:class:`Bin`\s are generated on the fly as needed.

  This implementation has two main benefits:

  * Decoupling of :cpp:class:`BinTable` resolution and memory requirements
  * Lookups in constant or linear time complexity with performance independent of resolution.

  **Constructors**

  .. cpp:function:: BinTable() = default;
  .. cpp:function:: BinTable(Reference chroms, std::uint32_t bin_size, std::size_t bin_offset = 0);
  .. cpp:function:: template <typename ChromIt> BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size, std::size_t bin_offset = 0);
  .. cpp:function:: template <typename ChromNameIt, typename ChromSizeIt> BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name, ChromSizeIt first_chrom_size, std::uint32_t bin_size, std::size_t bin_offset = 0);

  **Operators**

  .. cpp:function:: [[nodiscard]] bool operator==(const BinTable &other) const;
  .. cpp:function:: [[nodiscard]] bool operator!=(const BinTable &other) const;

  **Accessors**

  .. cpp:function:: [[nodiscard]] std::size_t size() const noexcept;
  .. cpp:function:: [[nodiscard]] bool empty() const noexcept;
  .. cpp:function:: [[nodiscard]] std::size_t num_chromosomes() const;
  .. cpp:function:: [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr const Reference &chromosomes() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr const std::vector<std::uint64_t> &num_bin_prefix_sum() const noexcept;

  **Iteration**

  .. cpp:function:: [[nodiscard]] auto begin() const -> iterator;
  .. cpp:function:: [[nodiscard]] auto end() const -> iterator;
  .. cpp:function:: [[nodiscard]] auto cbegin() const -> iterator;
  .. cpp:function:: [[nodiscard]] auto cend() const -> iterator;

  **Slicing**

  .. cpp:function:: [[nodiscard]] BinTable subset(const Chromosome &chrom) const;
  .. cpp:function:: [[nodiscard]] BinTable subset(std::string_view chrom_name) const;
  .. cpp:function:: [[nodiscard]] BinTable subset(std::uint32_t chrom_id) const;

  **Lookup**

  .. cpp:function:: [[nodiscard]] auto find_overlap(const GenomicInterval &query) const -> std::pair<BinTable::iterator, BinTable::iterator>;
  .. cpp:function:: [[nodiscard]] auto find_overlap(const Chromosome &chrom, std::uint32_t start, std::uint32_t end) const -> std::pair<BinTable::iterator, BinTable::iterator>;
  .. cpp:function:: [[nodiscard]] auto find_overlap(std::string_view chrom_name, std::uint32_t start, std::uint32_t end) const -> std::pair<BinTable::iterator, BinTable::iterator>;
  .. cpp:function:: [[nodiscard]] auto find_overlap(std::uint32_t chrom_id, std::uint32_t start, std::uint32_t end) const -> std::pair<BinTable::iterator, BinTable::iterator>;
  .. cpp:function:: [[nodiscard]] std::pair<Bin, Bin> at(const GenomicInterval &gi) const;
  .. cpp:function:: [[nodiscard]] std::pair<std::uint64_t, std::uint64_t> map_to_bin_ids(const GenomicInterval &gi) const;

  Query bins by genomic interval.

  .. cpp:function:: [[nodiscard]] Bin at(std::uint64_t bin_id) const;
  .. cpp:function:: [[nodiscard]] Bin at(const Chromosome &chrom, std::uint32_t pos = 0) const;
  .. cpp:function:: [[nodiscard]] Bin at(std::string_view chrom_name, std::uint32_t pos = 0) const;
  .. cpp:function:: [[nodiscard]] Bin at(std::uint32_t chrom_id, std::uint32_t pos) const;
  .. cpp:function:: [[nodiscard]] Bin at_hint(std::uint64_t bin_id, const Chromosome &chrom) const;

  Query by bin identifier.

  .. cpp:function:: [[nodiscard]] std::uint64_t map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const;
  .. cpp:function:: [[nodiscard]] std::uint64_t map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const;
  .. cpp:function:: [[nodiscard]] std::uint64_t map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const;

  Query by genomic coordinates

  **Others**

  .. cpp:function:: [[nodiscard]] BinTableConcrete concretize() const;

Pixels
------

.. cpp:namespace:: hictk

.. cpp:class:: template <typename N> ThinPixel

  Struct to model a genomic pixel using as little memory as possible.

  **Member variables**

  .. cpp:member:: static constexpr auto null_id = std::numeric_limits<std::uint64_t>::max();
  .. cpp:member:: std::uint64_t bin1_id{null_id};
  .. cpp:member:: std::uint64_t bin2_id{null_id};
  .. cpp:member:: N count{};

  **Factory methods**

  .. cpp:function:: static auto from_coo(std::string_view line) -> ThinPixel;
  .. cpp:function:: static auto from_coo(const BinTable &bins, std::string_view line) -> ThinPixel;

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator==(const ThinPixel &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const ThinPixel &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<(const ThinPixel &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<=(const ThinPixel &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>(const ThinPixel &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>=(const ThinPixel &other) const noexcept;


.. cpp:class:: PixelCoordinates;

  Struct to model 2D genomic coordinates using a pair of :cpp:class:`Bin`\s.

  **Member variables**

  .. cpp:member:: Bin bin1
  .. cpp:member:: Bin bin2

  **Constructors**

  .. cpp:function:: PixelCoordinates() = default;
  .. cpp:function:: PixelCoordinates(Bin bin1_, Bin bin2_) noexcept;
  .. cpp:function:: explicit PixelCoordinates(std::pair<Bin, Bin> bins) noexcept;
  .. cpp:function:: explicit PixelCoordinates(Bin bin) noexcept;

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator==(const PixelCoordinates &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const PixelCoordinates &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<(const PixelCoordinates &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<=(const PixelCoordinates &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>(const PixelCoordinates &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>=(const PixelCoordinates &other) const noexcept;

  **Accessors**

  .. cpp:function:: [[nodiscard]] bool is_intra() const noexcept;


.. cpp:class:: template <typename N> Pixel

  Struct to model genomic pixels as interaction counts associated to a pair of genomic :cpp:class:`Bin`\s.

  The main difference between :cpp:class:`ThinPixel` and :cpp:class:`Pixel` objects, is that the latter posesses all the knowledge required to map interactions to genomic coordinates, not just bin IDs.

  **Member variables**

  .. cpp:member:: PixelCoordinates coords{};
  .. cpp:member:: N count{};

  **Constructors**

  .. cpp:function:: Pixel() = default;
  .. cpp:function:: explicit Pixel(Bin bin, N count_ = 0) noexcept;
  .. cpp:function:: Pixel(Bin bin1_, Bin bin2_, N count_ = 0) noexcept;
  .. cpp:function:: explicit Pixel(PixelCoordinates coords_, N count_ = 0) noexcept;
  .. cpp:function:: Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end, N count_ = 0) noexcept;
  .. cpp:function:: Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1, const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2, N count_ = 0) noexcept;
  .. cpp:function:: Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_ = 0);
  .. cpp:function:: Pixel(const BinTable &bins, std::uint64_t bin_id, N count_ = 0);
  .. cpp:function:: Pixel(const BinTable &bins, const ThinPixel<N> &p);


  **Factory methods**

  .. cpp:function:: static auto from_coo(const BinTable &bins, std::string_view line) -> Pixel;
  .. cpp:function:: static auto from_bg2(const BinTable &bins, std::string_view line) -> Pixel;
  .. cpp:function:: static auto from_validpair(const BinTable &bins, std::string_view line) -> Pixel;
  .. cpp:function:: static auto from_4dn_pairs(const BinTable &bins, std::string_view line) -> Pixel;

  **Operators**

  .. cpp:function:: [[nodiscard]] explicit operator bool() const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator==(const Pixel<N> &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator!=(const Pixel<N> &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<(const Pixel<N> &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator<=(const Pixel<N> &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>(const Pixel<N> &other) const noexcept;
  .. cpp:function:: [[nodiscard]] bool operator>=(const Pixel<N> &other) const noexcept;

  **Conversion**

  .. cpp:function:: [[nodiscard]] ThinPixel<N> to_thin() const noexcept;
