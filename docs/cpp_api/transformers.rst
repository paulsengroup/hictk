..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

.. cpp:namespace:: hictk

Pixel transformers
==================

The transformer library provides a set of common algorithms used to manipulate streams of pixels.
Classes in defined in this library take a pair of pixel iterators or :cpp:class:`PixelSelector`\s directly and transform and/or aggregate them in different ways.

.. cpp:namespace:: hictk::transformers

Common types
------------

.. cpp:enum-class:: QuerySpan

  .. cpp:enumerator:: lower_triangle
  .. cpp:enumerator:: upper_triangle
  .. cpp:enumerator:: full

Coarsening pixels
-----------------

.. cpp:class:: template <typename PixelIt> CoarsenPixels

  Class used to coarsen pixels read from a pair of pixel iterators.
  The underlying sequence of pixels are expected to be sorted by their genomic coordinates.
  Coarsening is performed in a streaming fashion, minimizing the number of pixels that are kept into memory at any given time.

  .. cpp:function:: CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel,  std::shared_ptr<const BinTable> source_bins, std::size_t factor);

   Constructor for :cpp:class:`CoarsenPixels` class.
   ``first_pixel`` and ``last_pixels`` should be a pair of iterators pointing to the stream of pixels to be coarsened.
   ``source_bins`` is a shared pointer to the bin table to which ``first_pixel`` and ``last_pixel`` refer to.
   ``factor`` should be an integer value greater than 1, and is used to determine the properties of the ``target_bins`` :cpp:class:`BinTable` used for coarsening.

  **Accessors**

  .. cpp:function:: [[nodiscard]] const BinTable &src_bins() const noexcept;
  .. cpp:function:: [[nodiscard]] const BinTable &dest_bins() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> src_bins_ptr() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> dest_bins_ptr() const noexcept;

  :cpp:class:`BinTable` accessors.

  **Iteration**

  .. cpp:function:: [[nodiscard]] begin() const -> iterator;
  .. cpp:function:: [[nodiscard]] end() const -> iterator;
  .. cpp:function:: [[nodiscard]] cbegin() const -> iterator;
  .. cpp:function:: [[nodiscard]] cend() const -> iterator;

  Return an `InputIterator <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ to traverse the coarsened pixels.

  **Others**

  .. cpp:function:: [[nodiscard]] auto read_all() const -> std::vector<ThinPixel<N>>;

Selecting pixels overlapping with a band around the matrix diagonal
-------------------------------------------------------------------

.. cpp:class:: template <typename PixelIt> DiagonalBand

  Class used to select pixels overlapping with a band around the matrix diagonal.

  .. cpp:function:: DiagonalBand(PixelIt first_pixel, PixelIt last_pixel, std::uint64_t num_bins);

   Constructor for :cpp:class:`DiagonalBand` class.
   ``first_pixel`` and ``last_pixels`` should be a pair of iterators pointing to the stream of pixels to be processed.
   ``num_bins`` should correspond to the width of the band around the matrix diagonal.
   As all filtering operations are performed based on bin IDs, this transformer is unsuitable for processing pixels
   originating from files with bin tables of non-uniform size.

  **Iteration**

  .. cpp:function:: [[nodiscard]] begin() const -> iterator;
  .. cpp:function:: [[nodiscard]] end() const -> iterator;
  .. cpp:function:: [[nodiscard]] cbegin() const -> iterator;
  .. cpp:function:: [[nodiscard]] cend() const -> iterator;

  Return an `InputIterator <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ to traverse the pixels after filtering.

  **Others***

  .. cpp:function:: [[nodiscard]] auto read_all() const -> std::vector<Pixel<N>>;

Transforming COO pixels to BG2 pixels
-------------------------------------

.. cpp:class:: template <typename PixelIt> JoinGenomicCoords

  Class used to join genomic coordinates onto COO pixels, effectively transforming :cpp:class:`ThinPixel`\s into :cpp:class:`Pixel`\s.

  .. cpp:function:: JoinGenomicCoords(PixelIt first_pixel, PixelIt last_pixel,  std::shared_ptr<const BinTable> bins);

   Constructor for :cpp:class:`JoinGenomicCoords` class.
   ``first_pixel`` and ``last_pixels`` should be a pair of iterators pointing to the stream of pixels to be processed.
   ``bins`` is a shared pointer to the bin table to which ``first_pixel`` and ``last_pixel`` refer to.

  **Iteration**

  .. cpp:function:: [[nodiscard]] begin() const -> iterator;
  .. cpp:function:: [[nodiscard]] end() const -> iterator;
  .. cpp:function:: [[nodiscard]] cbegin() const -> iterator;
  .. cpp:function:: [[nodiscard]] cend() const -> iterator;

  Return an `InputIterator <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ to traverse the :cpp:class:`Pixel`\s.

  **Others***

  .. cpp:function:: [[nodiscard]] auto read_all() const -> std::vector<Pixel<N>>;


Merging streams of pre-sorted pixels
------------------------------------

.. cpp:class:: template <typename PixelIt> PixelMerger

  Class used to merge streams of pre-sorted pixels, yielding a sequence of unique pixels sorted by their genomic coordinates.
  Merging is performed in a streaming fashion, minimizing the number of pixels that are kept into memory at any given time.

  Duplicate pixels are aggregated by summing their corresponding interactions.
  Pixel merging also affects duplicate pixels coming from the same stream.

  .. cpp:function:: PixelMerger(std::vector<PixelIt> head, std::vector<PixelIt> tail);
  .. cpp:function:: template <typename ItOfPixelIt> PixelMerger(ItOfPixelIt first_head, ItOfPixelIt last_head, ItOfPixelIt first_tail);

  Constructors taking either two vectors of `InputIterators <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ or pairs of iterators to `InputIterators <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_.

  The ``head`` and ``tail`` vectors should contain the iterators pointing to the beginning and end of :cpp:class:`ThinPixel` streams, respectively.

  **Iteration**

  .. cpp:function:: [[nodiscard]] auto begin() const -> iterator;
  .. cpp:function:: [[nodiscard]] auto end() const noexcept -> iterator;

  Return an `InputIterator <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ to traverse the stream :cpp:class:`ThinPixel`\s after merging.

  **Others**

  .. cpp:function:: [[nodiscard]] auto read_all() const -> std::vector<PixelT>;


Computing common statistics
---------------------------

.. cpp:function:: template <typename PixelIt> [[nodiscard]] double avg(PixelIt first, PixelIt last);
.. cpp:function:: template <typename PixelIt, typename N> [[nodiscard]] N max(PixelIt first, PixelIt last);
.. cpp:function:: template <typename PixelIt> [[nodiscard]] std::size_t nnz(PixelIt first, PixelIt last);
.. cpp:function:: template <typename PixelIt, typename N> [[nodiscard]] N sum(PixelIt first, PixelIt last);


Converting streams of pixels to Arrow Tables
--------------------------------------------

.. cpp:enum-class:: DataFrameFormat

  .. cpp:enumerator:: COO
  .. cpp:enumerator:: BG2

.. cpp:class:: template <typename PixelIt> ToDataFrame

  .. cpp:function:: ToDataFrame(PixelIt first_pixel, PixelIt last_pixel, DataFrameFormat format = DataFrameFormat::COO, std::shared_ptr<const BinTable> bins = nullptr, QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false, bool mirror_pixels = true, std::size_t chunk_size = 256'000, std::optional<std::uint64_t> diagonal_band_width = {});
  .. cpp:function:: ToDataFrame(PixelIt first_pixel, PixelIt last_pixel, std::optional<PixelCoordinates> coord1_, std::optional<PixelCoordinates> coord2_ = {}, DataFrameFormat format = DataFrameFormat::COO, std::shared_ptr<const BinTable> bins = nullptr, QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false, bool mirror_pixels = true, std::size_t chunk_size = 256'000, std::optional<std::uint64_t> diagonal_band_width = {});
  .. cpp:function:: template <typename PixelSelector> ToDataFrame(const PixelSelector& sel, PixelIt it, DataFrameFormat format = DataFrameFormat::COO, std::shared_ptr<const BinTable> bins = nullptr, QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false, std::size_t chunk_size = 256'000, std::optional<std::uint64_t> diagonal_band_width = {});

  Construct an instance of a :cpp:class:`ToDataFrame` converter given a stream of pixels delimited by ``first_pixel`` and ``last_pixel``, a DataFrame ``format`` and a :cpp:class:`BinTable`.
  The underlying sequence of pixels are expected to be sorted by their genomic coordinates.

  The optional argument ``span`` determines whether the resulting :cpp:class:`arrow::Table` should contain interactions spanning the upper/lower-triangle or all interactions (regardless of whether they are located above or below the genome-wide matrix diagonal).
  It should be noted that queries spanning the the full-matrix or the lower-triangle are always more expensive because they involve an additional step where pixels are sorted by their genomic coordinates.

  When provided, the ``diagonal_band_width`` argument has the same semantics as the ``num_bins`` argument from the :cpp:class:`DiagonalBand` constructor.

  When fetching interactions with ``span=full`` from the cis portion of interaction maps using one of the overloads taking a pair of pixel iterators, users should provide a pair of genomic coordinates to ensure that, if necessary, interactions are correctly mirrored.

  .. cpp:function:: [[nodiscard]] std::shared_ptr<arrow::Table> operator()();

  Convert the stream of pixels into an :cpp:class:`arrow::Table`.


Converting streams of pixels to Eigen Dense Matrices
----------------------------------------------------

.. cpp:class:: template <typename N, typename PixelSelector> ToDenseMatrix

  .. cpp:type:: MatrixT = Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  .. cpp:function:: ToDenseMatrix(PixelSelector sel, N n, QuerySpan span = QuerySpan::full, std::optional<std::uint64_t> diagonal_band_width = {});
  .. cpp:function:: ToDenseMatrix(std::shared_ptr<const PixelSelector> sel, N n, QuerySpan span = QuerySpan::full, std::optional<std::uint64_t> diagonal_band_width = {});

  Construct an instance of a :cpp:class:`ToDenseMatrix` converter given a :cpp:class:`PixelSelector` object and a count type ``n``.

  The optional argument ``span`` determines whether the resulting matrix should contain interactions spanning the upper/lower-triangle or all interactions (regardless of whether they are located above or below the genome-wide matrix diagonal).
  Note that attempting to fetch trans-interactions with ``span=QuerySpan::lower_triangle`` will result in an exception being thrown.
  If you need to fetch trans-interactions from the lower-triangle, consider exchanging the range arguments used to fetch interactions, then transpose the resulting matrix.

  When provided, the ``diagonal_band_width`` argument has the same semantics as the ``num_bins`` argument from the :cpp:class:`DiagonalBand` constructor.

  .. cpp:function:: [[nodiscard]] auto operator()() -> MatrixT;

  Convert the stream of pixels into an :cpp:type:`MatrixT`.

Converting streams of pixels to Eigen Sparse Matrices
-----------------------------------------------------

.. cpp:class:: template <typename N, typename PixelSelector> ToSparseMatrix

  .. cpp:type:: MatrixT = Eigen::SparseMatrix<N, Eigen::RowMajor>;

  .. cpp:function:: ToSparseMatrix(PixelSelector sel, N n, QuerySpan span = QuerySpan::upper_triangle, bool minimize_memory_usage = false, std::optional<std::uint64_t> diagonal_band_width = {});
  .. cpp:function:: ToSparseMatrix(std::shared_ptr<const PixelSelector> sel, N n, QuerySpan span = QuerySpan::full, bool minimize_memory_usage = false, std::optional<std::uint64_t> diagonal_band_width = {});

  Construct an instance of a :cpp:class:`ToSparseMatrix` converter given a :cpp:class:`PixelSelector` object and a count type ``n``.

  The optional argument ``span`` determines whether the resulting matrix should contain interactions spanning the upper/lower-triangle or all interactions (regardless of whether they are located above or below the genome-wide matrix diagonal).
  Note that attempting to fetch trans-interactions with ``span=QuerySpan::lower_triangle`` will result in an exception being thrown.
  If you need to fetch trans-interactions from the lower-triangle, consider exchanging the range arguments used to fetch interactions, then transpose the resulting matrix.

  When ``minimize_memory_usage=true``, hictk will minimize memory usage by doing two passes over the queried pixels: one to calculate the exact number of entries to allocate for each row in the matrix, and the second pass to fill values in the matrix. This is usually slower than the default strategy, which traverses the data only once (but may overall require more memory than what is strictly needed). It should be noted that matrices are always compressed before being returned. Thus, the memory footprint of the matrices returned by :cpp:func:`ToSparseMatrix::operator()()` will be the same regardless of the fill strategy.

  When provided, the ``diagonal_band_width`` argument has the same semantics as the ``num_bins`` argument from the :cpp:class:`DiagonalBand` constructor.

  .. cpp:function:: [[nodiscard]] auto operator()() -> MatrixT;

  Convert the stream of pixels into an :cpp:type:`MatrixT`.
