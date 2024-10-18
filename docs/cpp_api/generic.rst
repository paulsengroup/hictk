..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

.. cpp:namespace:: hictk

Generic API
===========

hictk generic API allows users to transparently operate on .hic .cool files.
There is virtually no runtime overhead when using the :cpp:class:`File` and :cpp:class:`PixelSelector` classes. However iterating over :cpp:class:`Pixel`\s using this API is slightly slower than using the format-specific APIs.

Refer to examples in the :doc:`../quickstart_api` section for how to use the generic API without incurring into any overhead when iterating over :cpp:class:`Pixel`\s overlapping queries.

Common
------

.. cpp:namespace:: hictk

.. cpp:enum-class:: QUERY_TYPE

  .. cpp:enumerator:: BED
  .. cpp:enumerator:: UCSC


File handle
-----------

.. cpp:namespace:: hictk

.. cpp:class:: File

  This class implements a generic file handle capable of transparently operating on .cool and .hic files.

  **Constructors**

  .. cpp:function:: File(cooler::File clr);
  .. cpp:function:: File(hic::File hf);
  .. cpp:function:: File(std::string uri, std::uint32_t resolution = 0, hic::MatrixType type = hic::MatrixType::observed, hic::MatrixUnit unit = hic::MatrixUnit::BP);

   Constructors for :cpp:class:`File` class.
   ``resolution`` is a mandatory argument when opening .hic files.
   Matrix ``type`` and ``unit`` are ignored when operating on .cool files.

  **Accessors**

  .. cpp:function:: [[nodiscard]] std::string uri() const;

  Returns the URI of the open file. Always returns the file path when file is .hic.

  .. cpp:function:: [[nodiscard]] std::string path() const;

  Returns the path to the open file.

  .. cpp:function:: [[nodiscard]] constexpr bool is_hic() const noexcept;
  .. cpp:function:: [[nodiscard]] constexpr bool is_cooler() const noexcept;

  Test whether the open file is in .hic or .cool format.

  .. cpp:function:: [[nodiscard]] auto chromosomes() const -> const Reference &;
  .. cpp:function:: [[nodiscard]] auto bins() const -> const BinTable &;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const;

  Accessors to the chromosomes and bin table of the open file.

  .. cpp:function:: [[nodiscard]] std::uint32_t resolution() const;
  .. cpp:function:: [[nodiscard]] std::uint64_t nbins() const;
  .. cpp:function:: [[nodiscard]] std::uint64_t nchroms(bool include_ALL = false) const;

  Accessors for common attributes.
  Calling any of these accessors does not involve any computation.

  .. cpp:function:: [[nodiscard]] bool has_normalization(std::string_view normalization) const;
  .. cpp:function:: [[nodiscard]] std::vector<balancing::Method> avail_normalizations() const;
  .. cpp:function:: [[nodiscard]] const balancing::Weights &normalization(std::string_view normalization_) const;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(std::string_view normalization_) const;

  Accessors for normalization methods/vectors.

  **Fetch methods (1D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(const balancing::Method &normalization = balancing::Method::NONE()) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range, const balancing::Method &normalization = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start, std::uint32_t end, const balancing::Method &normalization = balancing::Method::NONE()) const;

   Return a :cpp:class:`PixelSelector` object that can be used to fetch pixels overlapping 1D (symmetric) queries.

   **Example usage:**

   .. code-block:: cpp

      hictk::File f{"myfile.hic", 1'000};

      // Fetch all pixels
      const auto sel1 = f.fetch();

      // Fetch all pixels (normalized with VC);
      const auto sel2 = f.fetch(balancing::Method::VC());

      // Fetch pixels overlapping chr1
      const auto sel3 = f.fetch("chr1");

      // Fetch pixels overlapping a region of interest
      const auto sel4 = f.fetch("chr1:10,000,000-20,000,000");
      const auto sel5 = f.fetch("chr1", 10'000'000, 20'000'000");

      // Fetch pixels using a BED query
      const auto sel6 = f.fetch("chr1\t10000000\t20000000",
                                balancing::Method::NONE(),
                                QUERY_TYPE::BED);

  **Fetch methods (2D queries)**

  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2, const balancing::Method &normalization = balancing::Method::NONE(), QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  .. cpp:function:: [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1, std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2, const balancing::Method &normalization = balancing::Method::NONE()) const;

  Return a :cpp:class:`PixelSelector` object that can be used to fetch pixels overlapping 2D (asymmetric) queries.

   **Example usage:**

   .. code-block:: cpp

      hictk::File f{"myfile.hic", 1'000};

      // Fetch pixels overlapping chr1:chr2
      const auto sel1 = f.fetch("chr1", "chr2");

      // Fetch pixels overlapping a region of interest
      const auto sel2 = f.fetch("chr1:10,000,000-20,000,000",
                                "chr2:10,000,000-20,000,000");
      const auto sel3 = f.fetch("chr1", 10'000'000, 20'000'000,
                                "chr2", 10'000'000, 20'000'000);


  **Advanced**

  .. cpp:function:: template <typename FileT> [[nodiscard]] constexpr const FileT &get() const;
  .. cpp:function:: template <typename FileT> [[nodiscard]] constexpr FileT &get();
  .. cpp:function:: [[nodiscard]] constexpr auto get() const noexcept -> const FileVar &;
  .. cpp:function:: [[nodiscard]] constexpr auto get() noexcept -> FileVar &;

   Methods to get the underlying :cpp:class:`hic::File` or :cpp:class:`cooler::File` file handle or a :cpp:class:`std::variant` of thereof.

   **Example usage:**

   .. code-block:: cpp

      hictk::File f{"myfile.hic", 1'000};

      assert(f.get<hic::File>().path() == "myfile.hic");
      assert(f.get<cooler::File>().path() == "myfile.hic");  // Throws an exception

      const auto fvar = f.get();
      std::visit([](const auto& f) {
        assert(f.path() == "myfile.hic");
      }, fvar);

Pixel selector
--------------

.. cpp:namespace:: hictk

.. cpp:class:: PixelSelector

  This class implements a generic, lightweight pixel selector object.

  :cpp:class:`PixelSelector` objects are constructed and returned by :cpp:func:`File::fetch` methods.
  Users are **not** supposed to construct :cpp:class:`PixelSelector` objects themselves.

  **Iteration**

  .. cpp:function:: template <typename N> [[nodiscard]] auto begin(bool sorted = true) const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto end() const -> iterator<N>;

  .. cpp:function:: template <typename N> [[nodiscard]] auto cbegin(bool sorted = true) const -> iterator<N>;
  .. cpp:function:: template <typename N> [[nodiscard]] auto cend() const -> iterator<N>;

  Return an `InputIterator <https://en.cppreference.com/w/cpp/named_req/InputIterator>`_ to traverse pixels
  overlapping the genomic coordinates used to create the :cpp:class:`PixelSelector`.

  Specifying ``sorted = false`` will improve throughput for queries over .hic files.

  When operating on .cool files, pixels are always returned sorted by genomic coordinates.

   **Example usage:**

   .. code-block:: cpp

      hictk::File f{"myfile.hic", 1'000};
      const auto sel = f.fetch();

      std::for_each(sel.begin<std::int32_t>(), sel.end<std::int32_t>(),
                    [&](const auto& pixel) { fmt::print("{}\n", pixel); });

      // STDOUT
      // 0  0 12
      // 0  2 7
      // 0  4 1
      // ...

  **Fetch at once**

  .. cpp:function:: template <typename N> [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  Read and return all :cpp:class:`Pixel`\s at once using a :cpp:class:`std::vector`.

  **Accessors**

  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  .. cpp:function:: [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  Return the genomic coordinates used to construct the :cpp:class:`PixelSelector`.

  .. cpp:function:: [[nodiscard]] const BinTable &bins() const noexcept;
  .. cpp:function:: [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;

  Return the :cpp:class:`BinTable` used to map :cpp:class:`Pixel`\s to genomic :cpp:class:`Bin`\s.

  .. cpp:function:: [[nodiscard]] const balancing::Weights &weights() const noexcept;

  Return the balancing weights associated with the :cpp:class:`PixelSelector` instance.

  **Advanced**

  .. cpp:function:: template <typename PixelSelectorT> [[nodiscard]] constexpr const PixelSelectorT &get() const;
  .. cpp:function:: template <typename PixelSelectorT> [[nodiscard]] constexpr PixelSelectorT &get();
  .. cpp:function:: [[nodiscard]] constexpr auto get() const noexcept -> const PixelSelectorVar &;
  .. cpp:function:: [[nodiscard]] constexpr auto get() noexcept -> PixelSelectorVar &;

   **Example usage:**

   .. code-block:: cpp

      hictk::File f{"myfile.hic", 1'000};

      const auto sel = f.fetch();

      assert(f.get<hic::PixelSelector>().matrix_type() == hic::MatrixType::observed");
      f.get<cooler::PixelSelector>();  // Throws an exception

      const auto selvar = sel.get();
      std::visit([](const auto& s) { assert(s.bins().resolution() == 1'000); }, selvar);
