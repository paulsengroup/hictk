..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Creating .cool and .hic files
#############################

hictk supports creating .cool and .hic files from text files in the following formats:

* `pairs (4DN-DCIC) <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md#example-pairs-file>`_
* `validPairs (nf-core/hic) <https://nf-co.re/hic/2.1.0/docs/output/#valid-pairs-detection-with-hic-pro>`_
* `bedGraph2 <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_
* `COO <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_

File requirements:

* ``dm6.chrom.sizes`` - `download <https://hgdownload.cse.ucsc.edu/goldenpath/dm6/bigZips/dm6.chrom.sizes>`__
* ``4DNFIKNWM36K.pairs.gz`` - `download <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/930ba072-05ac-4382-9a92-369517184ec7/4DNFIKNWM36K.pairs.gz>`__


Ingesting pairwise interactions into a 10kbp .cool file
-------------------------------------------------------

Loading interactions in pairs (4DN-DCIC) format into a .cool/hic file is straightforward:

.. code-block:: console

  user@dev:/tmp$ hictk load --format 4dn --bin-size 10kbp 4DNFIKNWM36K.pairs.gz 4DNFIKNWM36K.10000.cool

  [2024-09-26 16:51:28.059] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 16:51:28.068] [info]: begin loading pairwise interactions into a .cool file...
  [2024-09-26 16:51:28.137] [info]: writing chunk #1 to intermediate file "/tmp/hictk-tmp-XXXXQPdOSn/4DNFIKNWM36K.10000.cool.tmp"...
  [2024-09-26 16:51:45.281] [info]: done writing chunk #1 to tmp file "/tmp/hictk-tmp-XXXXQPdOSn/4DNFIKNWM36K.10000.cool.tmp".
  [2024-09-26 16:51:45.281] [info]: writing chunk #2 to intermediate file "/tmp/hictk-tmp-XXXXQPdOSn/4DNFIKNWM36K.10000.cool.tmp"...
  [2024-09-26 16:52:04.969] [info]: done writing chunk #2 to tmp file "/tmp/hictk-tmp-XXXXQPdOSn/4DNFIKNWM36K.10000.cool.tmp".
  [2024-09-26 16:52:04.970] [info]: merging 2 chunks into "4DNFIKNWM36K.10000.cool"...
  [2024-09-26 16:52:06.430] [info]: processing chr3L:1030000-1040000 chr3R:30240000-30250000 at 6882312 pixels/s...
  [2024-09-26 16:52:08.478] [info]: ingested 119208613 interactions (18122865 nnz) in 40.418916003s!

To ingest interactions in a .hic file, simply change the extension of the output file (or use the ``--output-fmt`` option).

hictk has native support for reading compressed interactions in the following formats: bzip2, lz4, lzo, gzip, xz, and zstd.

By default, the list of chromosomes is read from the file header.
The reference genome used to build the .cool or .hic file can be provided explicitly using the ``--chrom-sizes`` option.
Note that ``--chrom-sizes`` is a mandatory option when ingesting interactions in formats other than ``--format=4dn``.
In case the input file contains interactions mapping on chromosomes missing from the reference genome provided through ``--chrom-sizes``, the ``--drop-unknown-chroms`` flag can be used to instruct hictk to ignored said interactions.

When loading interactions using ``--format=pairs`` or ``--format=validPairs`` into a .cool file, tables of variable bins are supported.
To load interactions in to a .cool with a variable bin size provide the table of bins using the ``--bin-table`` option.

**Tips:**

* When creating large .cool/hic files, ``hictk`` needs to create potentially large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space.
* When loading interactions into .hic files, some of the steps can be run in parallel by increasing the number of processing threads using the ``--threads`` option.
* When loading pre-binned interactions into .cool file, if the interactions are already sorted by genomic coordinates, the ``--assume-sorted`` option can be used to load interactions at once, without using temporary files.
* Interaction loading performance can be improved by processing interactions in larger chunks. This can be controlled using the ``--chunk-size`` option. In fact, when ``--chunk-size`` is greater than the number of interactions to be loaded, .hic and .cool files can be created without the use of temporary files.


Merging multiple files
----------------------

Multiple .cool and .hic files using the same reference genome and resolution can be merged using ``hictk merge``:

.. code-block:: console

  # Merge multiple cooler files

  user@dev:/tmp$ hictk merge data/4DNFIZ1ZVXC8.mcool::/resolutions/10000 data/4DNFIZ1ZVXC8.mcool::/resolutions/10000 -o 4DNFIZ1ZVXC8.merged.10000.cool

  [2024-09-26 17:07:57.101] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 17:07:57.101] [info]: begin merging 2 files into one .cool file...
  [2024-09-26 17:07:58.978] [info]: processing chr3L:1030000-1040000 chr3R:29720000-29730000 at 5571031 pixels/s...
  [2024-09-26 17:08:01.224] [info]: DONE! Merging 2 files took 4.12s!
  [2024-09-26 17:08:01.224] [info]: data/4DNFIZ1ZVXC8.merged.10000.cool size: 19.64 MB

Merging .hic files as well as a mix of .hic and .cool files is also supported (as long as all files have the same resolution and reference genome).
When all input files contain data for multiple resolutions the ``--resolution`` option is mandatory.

**Tips:**

See the list of Tips for hictk load.
