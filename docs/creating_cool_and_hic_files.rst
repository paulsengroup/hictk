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

.. code-block:: console

  # Create a 10kbp .cool file using dm6 as reference
  user@dev:/tmp$ zcat 4DNFIKNWM36K.pairs.gz | hictk load --format 4dn --assembly dm6 --bin-size 10000 dm6.chrom.sizes 4DNFIKNWM36K.10000.cool

  [2024-01-23 15:15:00.520] [info]: Running hictk v0.0.6-45c36af-dirty
  [2024-01-23 15:15:00.531] [info]: writing chunk #1 to intermediate file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp"...
  [2024-01-23 15:15:23.762] [info]: done writing chunk #1 to tmp file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp".
  [2024-01-23 15:15:23.762] [info]: writing chunk #2 to intermediate file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp"...
  [2024-01-23 15:15:49.042] [info]: done writing chunk #2 to tmp file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp".
  [2024-01-23 15:15:49.042] [info]: writing chunk #3 to intermediate file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp"...
  [2024-01-23 15:15:49.834] [info]: done writing chunk #3 to tmp file "/tmp/4DNFIKNWM36K.10000.cool.tmp/4DNFIKNWM36K.10000.cool.tmp".
  [2024-01-23 15:15:49.836] [info]: merging 3 chunks into "4DNFIKNWM36K.10000.cool"...
  [2024-01-23 15:15:55.118] [info]: processing chr3L:15100000-15110000 chr3L:16230000-16240000 at 4789272 pixels/s...
  [2024-01-23 15:15:59.718] [info]: ingested 119208613 interactions (18122865 nnz) in 59.197723453s!

  # Create a 10kbp .hic file using dm6 as reference
  user@dev:/tmp$ zcat 4DNFIKNWM36K.pairs.gz | hictk load --format 4dn --assembly dm6 --bin-size 10000 dm6.chrom.sizes 4DNFIKNWM36K.10000.hic

  [2024-01-23 15:45:19.969] [info]: Running hictk v0.0.6-570037c-dirty
  [2024-01-23 15:45:42.439] [info]: preprocessing chunk #1 at 452919 pixels/s...
  [2024-01-23 15:46:09.182] [info]: preprocessing chunk #2 at 303750 pixels/s...
  [2024-01-23 15:46:11.184] [info]: writing header at offset 0
  [2024-01-23 15:46:11.184] [info]: begin writing interaction blocks to file "4DNFIKNWM36K.10000.hic"...
  [2024-01-23 15:46:11.184] [info]: [10000 bp] writing pixels for chr3R:chr3R matrix at offset 50632...
  [2024-01-23 15:46:13.295] [info]: [10000 bp] written 2264963 pixels for chr3R:chr3R matrix
  [2024-01-23 15:46:13.295] [info]: [10000 bp] writing pixels for chr3R:chr3L matrix at offset 4235718...
  [2024-01-23 15:46:14.611] [info]: [10000 bp] written 1610264 pixels for chr3R:chr3L matrix
  ...
  [2024-01-23 15:46:44.065] [info]: [10000 bp] initializing expected value vector
  [2024-01-23 15:46:50.531] [info]: [10000 bp] computing expected vector density
  [2024-01-23 15:46:51.157] [info]: writing 1 expected value vectors at offset 32065110...
  [2024-01-23 15:46:51.158] [info]: writing 0 normalized expected value vectors at offset 32078017...
  [2024-01-23 15:46:51.194] [info]: ingested 119208613 interactions (18122865 nnz) in 91.225341628s!

**Tips:**

* When creating large .hic files, ``hictk`` needs to create potentially large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space.


Merging multiple files
----------------------

Multiple .cool and .hic files using the same reference genome and resolution can be merged using ``hictk merge``:

.. code-block:: console

  # Merge multiple cooler files

  user@dev:/tmp$ hictk merge data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 -o 4DNFIZ1ZVXC8.merged.cool

  [2023-09-29 19:24:49.479] [info]: Running hictk v0.0.2
  [2023-09-29 19:24:49.479] [info]: begin merging 2 coolers...
  [2023-09-29 19:24:52.032] [info]: processing chr2R:11267000-11268000 chr4:1052000-1053000 at 3976143 pixels/s...
  [2023-09-29 19:24:55.157] [info]: processing chr3R:5812000-5813000 chr3R:23422000-23423000 at 3201024 pixels/s...
  [2023-09-29 19:24:57.992] [info]: DONE! Merging 2 coolers took 8.51s!
  [2023-09-29 19:24:57.992] [info]: 4DNFIZ1ZVXC8.merged.cool size: 36.23 MB

  # Merge multiple .hic files

  user@dev:/tmp$ hictk merge data/4DNFIZ1ZVXC8.hic9 data/4DNFIZ1ZVXC8.hic9 -o 4DNFIZ1ZVXC8.10000.merged.hic --resolution 10000

  [2024-01-23 15:49:23.248] [info]: Running hictk v0.0.6-570037c-dirty
  [2024-01-23 15:49:23.248] [info]: begin merging 2 .hic files...
  [2024-01-23 15:49:31.101] [info]: ingesting pixels at 1352814 pixels/s...
  [2024-01-23 15:49:37.777] [info]: writing header at offset 0
  [2024-01-23 15:49:37.777] [info]: begin writing interaction blocks to file "4DNFIZ1ZVXC8.10000.merged.hic"...
  [2024-01-23 15:49:37.777] [info]: [10000 bp] writing pixels for chr2L:chr2L matrix at offset 212...
  [2024-01-23 15:49:39.060] [info]: [10000 bp] written 1433133 pixels for chr2L:chr2L matrix
  [2024-01-23 15:49:39.060] [info]: [10000 bp] writing pixels for chr2L:chr2R matrix at offset 2619165...
  ...
  [2024-01-23 15:49:58.624] [info]: [10000 bp] initializing expected value vector
  [2024-01-23 15:50:05.276] [info]: [10000 bp] computing expected vector density
  [2024-01-23 15:50:05.276] [info]: writing 1 expected value vectors at offset 31936601...
  [2024-01-23 15:50:05.276] [info]: writing 0 normalized expected value vectors at offset 31949508...
  [2024-01-23 15:50:05.299] [info]: DONE! Merging 2 files took 42.05s!
  [2024-01-23 15:50:05.299] [info]: 4DNFIZ1ZVXC8.10000.merged.hic size: 31.95 MB

**Tips:**

* When merging many, large .hic files, ``hictk`` needs to create potentially large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space.
