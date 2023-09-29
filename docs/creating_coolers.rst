..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Creating .cool files
####################

hictk supports creating .cool files from text files in the following formats:

* `pairs (4DN-DCIC) <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md#example-pairs-file>`_
* `validPairs (nf-core/hic) <https://nf-co.re/hic/2.1.0/docs/output/#valid-pairs-detection-with-hic-pro>`_
* `bedGraph2 <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_
* `COO <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_

File requirements:

* ``dm6.chrom.sizes`` - `download <https://hgdownload.cse.ucsc.edu/goldenpath/dm6/bigZips/dm6.chrom.sizes>`_
* ``4DNFIKNWM36K.pairs.gz`` - `download <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/930ba072-05ac-4382-9a92-369517184ec7/4DNFIKNWM36K.pairs.gz>`_

.. code-block:: console

  # Create a 10kbp .cool file using dm6 as reference
  user@dev:/tmp$ zcat 4DNFIKNWM36K.pairs.gz | hictk load --format 4dn --assembly dms6 dm6.chrom.sizes 10000 4DNFIKNWM36K.1000.cool

  [2023-09-29 19:16:51.962] [info]: Running hictk v0.0.2
  [2023-09-29 19:16:51.962] [info]: begin loading un-sorted pairs...
  [2023-09-29 19:16:51.983] [info]: writing chunk #1 to intermediate file "4DNFIKNWM36K.1000.cool.tmp"...
  [2023-09-29 19:17:39.715] [info]: done writing chunk #1 to tmp file "4DNFIKNWM36K.1000.cool.tmp".
  [2023-09-29 19:17:39.715] [info]: writing chunk #2 to intermediate file "4DNFIKNWM36K.1000.cool.tmp"...
  [2023-09-29 19:17:39.719] [info]: done writing chunk #2 to tmp file "4DNFIKNWM36K.1000.cool.tmp".
  [2023-09-29 19:17:39.721] [info]: merging 2 chunks into "4DNFIKNWM36K.1000.cool"...
  [2023-09-29 19:17:41.716] [info]: processing chr3L:15100000-15110000 chr3L:16220000-16230000 at 5073567 pixels/s...


Merging multiple Cooler files
-----------------------------

Multiple .cool files using the same reference genome and resolution can be merged using ``hictk merge``:

.. code-block:: console

  user@dev:/tmp$ hictk merge data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 -o 4DNFIZ1ZVXC8.merged.cool

  [2023-09-29 19:24:49.479] [info]: Running hictk v0.0.2
  [2023-09-29 19:24:49.479] [info]: begin merging 2 coolers...
  [2023-09-29 19:24:52.032] [info]: processing chr2R:11267000-11268000 chr4:1052000-1053000 at 3976143 pixels/s...
  [2023-09-29 19:24:55.157] [info]: processing chr3R:5812000-5813000 chr3R:23422000-23423000 at 3201024 pixels/s...
  [2023-09-29 19:24:57.992] [info]: DONE! Merging 2 coolers took 8.51s!
  [2023-09-29 19:24:57.992] [info]: 4DNFIZ1ZVXC8.merged.cool size: 36.23 MB
