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

  @article {hictk,
	  author = {Roberto Rossini and Jonas Paulsen},
	  title = {hictk: blazing fast toolkit to work with .hic and .cool files},
	  elocation-id = {2023.11.26.568707},
	  year = {2023},
	  doi = {10.1101/2023.11.26.568707},
	  publisher = {Cold Spring Harbor Laboratory},
	  URL = {https://www.biorxiv.org/content/early/2023/11/27/2023.11.26.568707},
	  eprint = {https://www.biorxiv.org/content/early/2023/11/27/2023.11.26.568707.full.pdf},
	  journal = {bioRxiv}
  }

If you use ``hictk convert`` to convert .[m]cool files to .hic format you should also cite JuicerTools or HiCTools.


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
