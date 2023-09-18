..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Introduction
============

hictk is a blazing fast toolkit to work with .hic and .cool files.

hictk is capable of reading files in .cool, .mcool, .scool and .hic format (including hic v9) as well as writing .cool, .mcool and .scool files.

.. only:: not latex

   Documentation formats
   ---------------------

   You are reading the HTML version of the documentation. An alternative `PDF
   version <https://hictk.readthedocs.io/_/downloads/en/latest/pdf/>`__ is
   also available.

   Installation
   ------------

.. only:: latex

   .. rubric:: Installation

hictk is developed on Linux and tested on Linux, MacOS and Windows. CLI tools can be installed in several different ways. Refer to :doc:`Installation <./installation>` for more details.

hictk can be compiled on most UNIX-like systems (including many Linux distributions and MacOS) as well as Windows. See the :doc:`build instructions <./installation_src>` for more details.

Python bindings for hictk can be installed using pip or conda. See :doc:`here <./installation_py>` for more details.

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite nanobind in scientific
discourse:

.. code-block:: bibtex

    @misc{hictk,
       author = {Roberto Rossini},
       year = {2023},
       note = {https://github.com/paulsengroup/hictk},
       title = {hictk: blazing fast toolkit to work with .hic and .cool files}
    }


.. only:: not latex

   Table of contents
   -----------------

.. toctree::
   :caption: Installation
   :maxdepth: 1

   installation
   installation_src
   installation_py

.. toctree::
   :caption: Introduction
   :maxdepth: 1

   quickstart
   file_validation


.. toctree::
   :caption: CLI and API Reference
   :maxdepth: 2

   cli_reference
   cpp_api
   cpp_api/index
   Python API <https://github.com/paulsengroup/hictkpy>
