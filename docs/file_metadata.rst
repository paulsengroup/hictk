..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

File metadata
#############

``hictk metadata`` can be used to read metadata from any Cooler file (that is, .cool, .mcool, and .scool files) as well as .hic files.

Fetching metadata from single-resolution files
----------------------------------------------

The following shows how to use ``hictk metadata`` to access metadata information stored in single-resolution Cooler files (or URIs, like shown in the example below):

.. code-block:: console

  user@dev:/tmp$ hictk metadata data/4DNFIZ1ZVXC8.mcool::/resolutions/1000

  {
      "assembly": "dm6",
      "bin-size": 1000,
      "bin-type": "fixed",
      "creation-date": "2023-07-07T14:15:30.759988",
      "format": "HDF5::Cooler",
      "format-url": "https://github.com/4dn-dcic/hic2cool",
      "format-version": 3,
      "generated-by": "hic2cool-0.8.3",
      "nbins": 137572,
      "nchroms": 8,
      "nnz": 26591454,
      "storage-mode": "symmetric-upper"
  }

By default, ``hictk metadata`` outputs information in JSON format, however the output format can be changed using the ``--output-format`` CLI options (currently, ``json``, ``toml``, and ``yaml`` formats are supported).


Fetching metadata from multi-resolution files
---------------------------------------------

Next, we look how to fetch metadata from multi-resolution .hic and .mcool files.

.. code-block:: console

  user@dev:/tmp$ hictk metadata data/4DNFIZ1ZVXC8.hic9

  {
      "assembly": "dm6",
      "format": "HIC",
      "format-url": "https://github.com/aidenlab/hic-format",
      "format-version": 9,
      "hicFileScalingFactor": 1.0,
      "nchroms": 8,
      "resolutions": [
          1000,
          5000,
          10000,
          25000,
          50000,
          100000,
          250000,
          500000,
          1000000,
          2500000
      ],
      "software": "Juicer Tools Version 3.30.00"
  }

When dealing with multi-resolution and single-cell files, it is possible to also view metadata information of individual resolutions/cells by using the ``--recursive`` CLI flag:

.. code-block:: console

  user@dev:/tmp$ hictk metadata data/4DNFIZ1ZVXC8.mcool --recursive

  {
      "1000": {
          "assembly": "dm6",
          "bin-size": 1000,
          "bin-type": "fixed",
          "creation-date": "2023-07-07T14:15:30.759988",
          "format": "HDF5::Cooler",
          "format-url": "https://github.com/4dn-dcic/hic2cool",
          "format-version": 3,
          "generated-by": "hic2cool-0.8.3",
          "nbins": 137572,
          "nchroms": 8,
          "nnz": 26591454,
          "storage-mode": "symmetric-upper"
      },
      ...
      "bin-type": "fixed",
      "format": "HDF5::MCOOL",
      "format-version": 2,
      "resolutions": [
          1000,
          5000,
          10000,
          25000,
          50000,
          100000,
          250000,
          500000,
          1000000,
          2500000
      ]
  }
