..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Format conversion
#################

``hictk`` supports conversion between .hic and .[m]cool file formats (including .hic v9 files).

Converting from .hic to .[m]cool
--------------------------------

Converting from .hic to .cool or .mcool formats consists of the following operations

#. Fetch the list of available resolutions
#. For each resolution to be converted:

   a. Copy all raw interactions present in the .hic file
   b. Copy all normalization vectors

Interactions are copied using streams of data, so memory requirements remain quite modest even when converting very high resolutions.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.hic9 4DNFIZ1ZVXC8.mcool

  [2024-09-26 16:06:41.713] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 16:06:41.713] [info]: Converting data/4DNFIZ1ZVXC8.hic9 to 4DNFIZ1ZVXC8.mcool (hic -> mcool)...
  [2024-09-26 16:06:41.943] [info]: [1000] begin processing 1000bp matrix...
  [2024-09-26 16:06:44.117] [info]: [1000] processing chr2R:11267000-11268000 at 4604052 pixels/s (cache hit rate 0.00%)...
  [2024-09-26 16:06:46.026] [info]: [1000] processing chr3R:5812000-5813000 at 5238345 pixels/s (cache hit rate 0.10%)...
  [2024-09-26 16:06:47.842] [info]: [1000] processing SCALE normalization vector...
  [2024-09-26 16:06:47.873] [info]: [1000] processing VC normalization vector...
  [2024-09-26 16:06:47.907] [info]: [1000] processing VC_SQRT normalization vector...
  [2024-09-26 16:06:48.411] [info]: [1000] DONE! Processed 26682908 pixels across 8 chromosomes in 6.47s
  ...
  [2024-09-26 16:06:58.265] [info]: DONE! Processed 10 resolution(s) in 16.55s!
  [2024-09-26 16:06:58.265] [info]: data/4DNFIZ1ZVXC8.hic9 size: 133.68 MB
  [2024-09-26 16:06:58.265] [info]: 4DNFIZ1ZVXC8.mcool size: 99.86 MB


It is also possible to convert only a subset of available resolutions by specifying resolutions to be converted with the ``--resolutions`` option.

When specifying a single resolution, the resulting file will be in .cool format.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.hic9 4DNFIZ1ZVXC8.1000.cool --resolutions 1000

  [2024-09-26 16:08:09.827] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 16:08:09.827] [info]: Converting data/4DNFIZ1ZVXC8.hic9 to 4DNFIZ1ZVXC8.cool (hic -> cool)...
  [2024-09-26 16:08:10.043] [info]: [1000] begin processing 1000bp matrix...
  [2024-09-26 16:08:11.216] [info]: [1000] processing chr2R:11267000-11268000 at 8539710 pixels/s (cache hit rate 93.05%)...
  [2024-09-26 16:08:12.462] [info]: [1000] processing chr3R:5812000-5813000 at 8032129 pixels/s (cache hit rate 93.11%)...
  [2024-09-26 16:08:13.423] [info]: [1000] processing SCALE normalization vector...
  [2024-09-26 16:08:13.453] [info]: [1000] processing VC normalization vector...
  [2024-09-26 16:08:13.485] [info]: [1000] processing VC_SQRT normalization vector...
  [2024-09-26 16:08:13.968] [info]: [1000] DONE! Processed 26682908 pixels across 8 chromosomes in 3.92s
  [2024-09-26 16:08:13.968] [info]: DONE! Processed 1 resolution(s) in 4.14s!
  [2024-09-26 16:08:13.968] [info]: data/4DNFIZ1ZVXC8.hic9 size: 133.68 MB
  [2024-09-26 16:08:13.968] [info]: 4DNFIZ1ZVXC8.cool size: 36.69 MB


Converting from .[m]cool to .hic
--------------------------------

``hictk convert`` can also be used to convert .[m]cool files to .hic format.

The conversion steps are similar to those carried out to convert .hic to .[m]cool.
The main difference is that in this case hictk computes the raw and normalized expected values for each resolution.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.mcool 4DNFIZ1ZVXC8.hic

  [2024-09-26 16:10:58.066] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 16:10:58.066] [info]: Converting data/4DNFIZ1ZVXC8.mcool to 4DNFIZ1ZVXC8.hic (mcool -> hic)...
  [2024-09-26 16:11:02.124] [info]: ingesting pixels at 2472799 pixels/s...
  [2024-09-26 16:11:06.328] [info]: ingesting pixels at 2379253 pixels/s...
  [2024-09-26 16:11:13.161] [info]: ingesting pixels at 2479544 pixels/s...
  [2024-09-26 16:11:17.436] [info]: ingesting pixels at 2339729 pixels/s...
  [2024-09-26 16:11:24.176] [info]: ingesting pixels at 2472188 pixels/s...
  [2024-09-26 16:11:32.941] [info]: writing header at offset 0
  [2024-09-26 16:11:32.941] [info]: begin writing interaction blocks to file "4DNFIZ1ZVXC8.hic"...
  [2024-09-26 16:11:32.941] [info]: [1000 bp] writing pixels for chr2L:chr2L matrix at offset 249...
  [2024-09-26 16:11:35.129] [info]: [1000 bp] written 2676654 pixels for chr2L:chr2L matrix
  [2024-09-26 16:11:35.159] [info]: [5000 bp] writing pixels for chr2L:chr2L matrix at offset 4075891...
  [2024-09-26 16:11:37.035] [info]: [5000 bp] written 2676654 pixels for chr2L:chr2L matrix
  [2024-09-26 16:11:37.096] [info]: [10000 bp] writing pixels for chr2L:chr2L matrix at offset 8697885...
  [2024-09-26 16:11:38.094] [info]: [10000 bp] written 1433133 pixels for chr2L:chr2L matrix
  ...
  [2024-09-26 16:13:20.981] [info]: [2500000 bp] initializing expected value vector
  [2024-09-26 16:13:20.981] [info]: [2500000 bp] computing expected vector density
  [2024-09-26 16:13:20.982] [info]: [500000 bp] computing expected vector density
  [2024-09-26 16:13:20.982] [info]: writing 50 normalized expected value vectors at offset 135622984...
  [2024-09-26 16:13:20.983] [info]: writing 400 normalization vectors at offset 136510590...
  [2024-09-26 16:13:21.027] [info]: DONE! Processed 10 resolution(s) in 142.96s!
  [2024-09-26 16:13:21.027] [info]: data/4DNFIZ1ZVXC8.mcool size: 139.37 MB
  [2024-09-26 16:13:21.027] [info]: 4DNFIZ1ZVXC8.hic size: 140.32 MB

**Tips:**

* When converting large .[m]cool files to .hic, ``hictk`` may need to create large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space.
* When converting .[m]cool files to .hic certain conversion steps can be performed in parallel. To improve performance, please make sure to increase the number of processing threads with option ``--threads``.
