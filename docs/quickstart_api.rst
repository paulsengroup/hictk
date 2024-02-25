..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Quickstart (API)
################

The library component of hictk, libhictk, can be installed and configured in several ways.

Installing libhictk
===================

Installing using Conan
----------------------

To install libhictk using Conan, first create a conanfile.txt like the following:

.. code-block::

  [requires]
  hictk/0.0.9

  [generators]
  CMakeDeps

  [layout]
  cmake_layout

Next, install hictk as follows:

.. code-block:: sh

  conan install conanfile.txt --build=missing --output-folder=conan_deps

Folder ``conan_deps`` will contain all CMake module and config files required to include hictk in an application using CMake as build generator.

Finally, add ``find_package(hictk REQUIRED)`` to your ``CMakeLists.txt`` and pass the full path to folder ``conan_deps`` to CMake through the ``CMAKE_PREFIX_PATH`` variable:

.. code-block:: sh

  cmake -DCMAKE_PREFIX_PATH='/path/to/conan_deps' ... -B build/ -S .


For more options and details refer to hictk page on `ConanCenter <https://conan.io/center/recipes/hictk>`_.

Installing using CMake FetchContent
-----------------------------------

Before beginning, make sure all of hictk dependencies have been installed.
Refer to `conanfile.txt <https://github.com/paulsengroup/hictk/blob/main/conanfile.txt>`_ for an up-to-date list of hictk dependencies.

To install and configure hictk using `FetchContent <https://cmake.org/cmake/help/latest/module/FetchContent.html>`_, first write a ``CMakeLists.txt`` file like the following:

.. code-block:: cmake

  cmake_minimum_required(VERSION 3.25)
  cmake_policy(VERSION 3.25...3.27)

  project(myproject LANGUAGES C CXX)

  include(FetchContent)
  FetchContent_Declare(
    hictk
    GIT_REPOSITORY  "https://github.com/paulsengroup/hictk.git"
    GIT_TAG         v0.0.9
    SYSTEM)

  # Customize hictk build flags
  set(HICTK_ENABLE_TESTING OFF)
  set(HICTK_BUILD_EXAMPLES OFF)
  set(HICTK_BUILD_BENCHMARKS OFF)
  set(HICTK_BUILD_TOOLS OFF)
  set(HICTK_INSTALL OFF)

  FetchContent_MakeAvailable(hictk)

  add_executable(main main.cpp)
  target_link_libraries(main PRIVATE hictk::file)  # Add other targets as necessary

Include hictk source using CMake add_subdirectory
-------------------------------------------------

Simply add a copy of hictk source code to your source tree (e.g. under folder ``myproject/external/hictk``), then add ``add_subdirectory("external/hictk")`` to your ``CMakeLists.txt``.


A quick tour of libhictk
------------------------

libhictk is a C++17 header-only library that provides the building blocks required to build complex applications operating on .hic and .cool files.

libhictk public API is organized in 5 main sections:

.. cpp:namespace:: hictk

#. Classes :cpp:class:`cooler::File`, :cpp:class:`cooler::MultiResFile` and :cpp:class:`cooler::SingleCellFile`, which can be used to read and write .cool, .mcool and .scool files respectively.
#. Class :cpp:class:`hic::File` which can be used to read .hic files
#. Class :cpp:class:`File` which wraps :cpp:class:`cooler::File` and :cpp:class:`hic::File` and provides a uniform interface to read .cool and .hic files
#. Various other classes used e.g. to model tables of bins, reference genomes and much more
#. Classes and free-standing functions to perform common operations on files or pixel iterators, such as coarsening and balancing.

The quick tour showcases basic functionality of the generic :cpp:class:`File` class. For more detailed examples refer to hictk `examples <https://github.com/paulsengroup/hictk/tree/main/examples>`_ and :doc:`cpp_api/index`.

.. code-block:: cpp

  #include <algorithm>
  #include <cstdint>
  #include <hictk/file.hpp>
  #include <iostream>
  #include <string>

  int main() {
    // const std::string path = "interactions.cool";
    // const std::string path = "interactions.mcool::/resolutions/1000";
    const std::string path = "interactions.hic";
    const std::uint32_t resolution = 1000;

    const hictk::File f(path, resolution);

    const auto selector = f.fetch("chr1", "chr2");

    std::for_each(selector.template begin<std::int32_t>(), selector.template end<std::int32_t>(),
                  [](const hictk::ThinPixel<std::int32_t>& p) {
                    std::cout << p.bin1_id << "\t";
                    std::cout << p.bin2_id << "\t";
                    std::cout << p.count << "\n";
                  });
  }


It is often the case that we need access to more information than just bin IDs and counts.
Joining genomic coordinates to pixel counts can be done as follows:

.. code-block:: cpp

  #include <cstdint>
  #include <hictk/file.hpp>
  #include <hictk/transformers.hpp>
  #include <iostream>
  #include <string>

  int main() {
    const std::string path = "interactions.hic";
    const std::uint32_t resolution = 1000;

    const hictk::File f(path, resolution);

    const auto selector = f.fetch("chr1", "chr2");
    const hictk::transformers::JoinGenomicCoords jselector(
        selector.template begin<std::int32_t>(), selector.template end<std::int32_t>(), f.bins_ptr());

    for (const auto& p : jselector) {
      std::cout << p.coords.bin1.chrom().name() << "\t";
      std::cout << p.coords.bin1.start() << "\t";
      std::cout << p.coords.bin1.end() << "\t";
      std::cout << p.coords.bin2.chrom().name() << "\t";
      std::cout << p.coords.bin2.start() << "\t";
      std::cout << p.coords.bin2.end() << "\t";
      std::cout << p.count << "\n";
    }
  }


The above examples work just fine, however using iterators returned by generic :cpp:class:`PixelSelector` is suboptimal. These iterators are implemented using `std::variant <https://en.cppreference.com/w/cpp/utility/variant>`_ and require checking the type of the underlying ``PixelSelector`` every iteration. The overhead of this check is quite low but still noticeable.

We can avoid paying this overhead by using the format-specific file handles instead of the generic one, or by using `std::visit <https://en.cppreference.com/w/cpp/utility/variant/visit>`_ as follows:

.. code-block:: cpp

  #include <algorithm>
  #include <cstdint>
  #include <hictk/file.hpp>
  #include <iostream>
  #include <string>
  #include <variant>

  int main() {
    const std::string path = "interactions.hic";
    const std::uint32_t resolution = 1000;

    const hictk::File f(path, resolution);

    const auto selector = f.fetch("chr1", "chr2");

    // std::visit applies the lambda function provided as first argument
    // to the variant returned by selector.get().
    // In this way, the type held by the std::variant is checked once
    // and the underlying PixelSelector and iterators are used for all operations
    std::visit(
        [&](const auto& sel) {
          std::for_each(sel.template begin<std::int32_t>(), sel.template end<std::int32_t>(),
                        [](const hictk::ThinPixel<std::int32_t>& p) {
                          std::cout << p.bin1_id << "\t";
                          std::cout << p.bin2_id << "\t";
                          std::cout << p.count << "\n";
                        });
        },
        selector.get());
  }
