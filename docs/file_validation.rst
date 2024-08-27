..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

File validation
###############

Why is this needed?
-------------------

``hictk validate`` can detect several types of data corruption in .hic and .cool files, from simple file truncation due to e.g. failed downloads to subtle index corruption in .cool files.

.. _cooler-index-corruption-label:

Cooler index corruption
^^^^^^^^^^^^^^^^^^^^^^^

To make a long story short, older versions of cooler (including v0.8.3) had a bug in ``cooler zoomify`` that caused the generation of invalid file indexes. This results in duplicate pixels with different values being reported for the affected region.

Example:

.. csv-table:: Output of cooler dump for corrupted file `4DNFI9GMP2J8.mcool <https://data.4dnucleome.org/files-processed/4DNFI9GMP2J8/>`_
  :file: ./assets/corrupted_mcool_example.tsv
  :header-rows: 1
  :delim: tab

Unfortunately, this is not a rare issue, as the above bug currently affects most (possibly all) .mcool files released by 4DNucleome:

.. image:: assets/4dnucleome_bug_notice.avif

hictk validate
--------------

``hictk validate`` was initially developed to detect files affected by the above issue and was later extended to also validate .cool, .scool and .hic files.

Perform a quick check to detect truncated or otherwise invalid files:

.. code-block:: console

  # Validate a .hic file
  user@dev:/tmp$ hictk validate test/data/hic/4DNFIZ1ZVXC8.hic8
  ### SUCCESS: "test/data/hic/4DNFIZ1ZVXC8.hic8" is a valid .hic file.

  # Validate a .cool file
  user@dev:/tmp$ hictk validate test/data/integration_tests/4DNFIZ1ZVXC8.mcool
  uri="test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/2500000"
  is_hdf5=true
  unable_to_open_file=false
  file_was_properly_closed=true
  missing_or_invalid_format_attr=false
  missing_or_invalid_bin_type_attr=false
  missing_groups=[]
  is_valid_cooler=true
  index_is_valid=not_checked
  ### SUCCESS: "test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/2500000" is a valid Cooler.
  uri="test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/1000000"
  is_hdf5=true
  unable_to_open_file=false
  file_was_properly_closed=true
  missing_or_invalid_format_attr=false
  missing_or_invalid_bin_type_attr=false
  missing_groups=[]
  is_valid_cooler=true
  index_is_valid=not_checked
  ### SUCCESS: "test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/1000000" is a valid Cooler.
  ...
  uri="test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/1000"
  is_hdf5=true
  unable_to_open_file=false
  file_was_properly_closed=true
  missing_or_invalid_format_attr=false
  missing_or_invalid_bin_type_attr=false
  missing_groups=[]
  is_valid_cooler=true
  index_is_valid=not_checked
  ### SUCCESS: "test/data/integration_tests/4DNFIZ1ZVXC8.mcool::/resolutions/1000" is a valid Cooler.


The quick check will not detect Cooler files with corrupted index, as this requires the ``--validate-index`` option:

.. code-block:: console

  user@dev:/tmp$ hictk validate --validate-index 4DNFI9GMP2J8.mcool::/resolutions/1000000
  uri="4DNFI9GMP2J8.mcool::/resolutions/1000000"
  is_hdf5=true
  unable_to_open_file=false
  file_was_properly_closed=true
  missing_or_invalid_format_attr=false
  missing_or_invalid_bin_type_attr=false
  missing_groups=[]
  is_valid_cooler=true
  index_is_valid=false
  ### FAILURE: "4DNFI9GMP2J8.mcool::/resolutions/1000000" is not a valid Cooler.

Restoring corrupted .mcool files
--------------------------------

Luckily, the base resolution of .mcool files corrupted as described in :ref:`cooler-index-corruption-label` is still valid, and so corrupted resolutions can be regenerated from the base resolution.

File restoration is automated with ``hictk fix-mcool``:

.. code-block:: sh

  hictk fix-mcool 4DNFI9GMP2J8.mcool 4DNFI9GMP2J8.fixed.mcool

``hictk fix-mcool`` is basically a wrapper around ``hictk zoomify`` and ``hictk balance``.

When balancing, ``hictk fix-mcool`` will try to use the same parameters used to balance the original .mcool file. When this is not possible, ``hictk fix-mcool`` will fall back to the default parameters used by ``hictk balance``.
