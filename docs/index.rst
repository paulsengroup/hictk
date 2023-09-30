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
   version <https://hictk.readthedocs.io/_/downloads/en/latest/pdf/>`_ is
   also available.

   Installation
   ------------

.. only:: latex

   .. rubric:: Installation

hictk is developed on Linux and tested on Linux, MacOS and Windows. CLI tools can be installed in several different ways. Refer to :doc:`Installation <./installation>` for more details.

hictk can be compiled on most UNIX-like systems (including many Linux distributions and MacOS) as well as Windows. See the :doc:`build instructions <./installation_src>` for more details.

Python bindings for hictk can be installed using pip or conda. Refer to hictkpy `documentation <https://hictkpy.readthedocs.io/en/latest/installation.html>`_ for more details.

.. only:: not latex

   How to cite this project?
   -------------------------

.. only:: latex

   .. rubric:: How to cite this project?

Please use the following BibTeX template to cite hictk in scientific
discourse:

.. code-block:: bibtex

    @misc{hictk,
       author = {Roberto Rossini},
       year = {2023},
       note = {https://github.com/paulsengroup/hictk},
       title = {hictk: blazing fast toolkit to work with .hic and .cool files}
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
   creating_coolers
   creating_multires_coolers


.. toctree::
   :caption: CLI and API Reference
   :maxdepth: 2

   cli_reference
   cpp_api/index
   Python API <https://hictkpy.readthedocs.io/en/latest/hictkpy.html>
