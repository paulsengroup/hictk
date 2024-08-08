..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Introduction
============

hictk is a blazing fast toolkit to work with .hic and .cool files.

hictk is capable of reading and writing files in .cool, .mcool, .scool and .hic format (including hic v9).

.. only:: not latex

   Documentation formats
   ---------------------

   You are reading the HTML version of the documentation. An alternative `PDF
   version <https://hictk.readthedocs.io/_/downloads/en/latest/pdf/>`_ is
   also available.

   Installation
   ------------

.. only:: latex

   .. rubric:: Installation

hictk is developed on Linux and tested on Linux, MacOS and Windows. CLI tools can be installed in several different ways. Refer to :doc:`Installation <./installation>` for more details.

hictk can be compiled on most UNIX-like systems (including many Linux distributions and MacOS) as well as Windows. See the :doc:`build instructions <./installation_src>` for more details.

hictk can be used from languages other than C++ through the following bindings:

* Python bindings through `hictkpy <https://github.com/paulsengroup/hictkpy>`_
* R bindings through `hictkR <https://github.com/paulsengroup/hictkR>`_

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite hictk in scientific
discourse:

.. code-block:: bibtex

  @article{hictk,
      author = {Rossini, Roberto and Paulsen, Jonas},
      title = "{hictk: blazing fast toolkit to work with .hic and .cool files}",
      journal = {Bioinformatics},
      volume = {40},
      number = {7},
      pages = {btae408},
      year = {2024},
      month = {06},
      issn = {1367-4811},
      doi = {10.1093/bioinformatics/btae408},
      url = {https://doi.org/10.1093/bioinformatics/btae408},
      eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/7/btae408/58385157/btae408.pdf},
  }


.. only:: not latex

   Table of contents
   -----------------

.. toctree::
   :caption: Installation
   :maxdepth: 1

   installation
   installation_src

.. toctree::
   :caption: Introduction
   :maxdepth: 1

   quickstart_cli
   quickstart_api
   downloading_test_datasets
   file_validation
   format_conversion
   reading_interactions
   creating_cool_and_hic_files
   creating_multires_files
   balancing_matrices

.. toctree::
   :caption: Tutorials
   :maxdepth: 2

   tutorials/reordering_chromosomes
   tutorials/dump_interactions_to_cool_hic_file


.. toctree::
   :caption: CLI and API Reference
   :maxdepth: 2

   cli_reference
   cpp_api/index
   Python API <https://hictkpy.readthedocs.io/en/latest/index.html>
   R API <https://paulsengroup.github.io/hictkR/reference/index.html>
