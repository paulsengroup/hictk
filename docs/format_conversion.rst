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
   b. Copy all known normalization vectors (currently these are VC, VC_SQRT, KR, and SCALE)

Interactions are copied using streams of data, so memory requirements remain quite modest even when converting very high resolutions.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.hic9 4DNFIZ1ZVXC8.mcool

  [2023-09-29 17:12:08.983] [info]: Running hictk v0.0.2-f83f93e
  [2023-09-29 17:12:08.983] [info]: Converting data/4DNFIZ1ZVXC8.hic9 to 4DNFIZ1ZVXC8.mcool (hic -> mcool)...
  [2023-09-29 17:12:09.052] [info]: [1000] begin processing 1000bp matrix...
  [2023-09-29 17:12:12.212] [info]: [1000] processing chr2R:11267000-11268000 at 3167564 pixels/s (cache hit rate 0.00%)...
  [2023-09-29 17:12:15.346] [info]: [1000] processing chr3R:5672000-5673000 at 3190810 pixels/s (cache hit rate 0.00%)...
  [2023-09-29 17:12:18.204] [info]: [1000] processing SCALE normalization vector...
  [2023-09-29 17:12:18.241] [info]: [1000] processing VC normalization vector...
  [2023-09-29 17:12:18.285] [info]: [1000] processing VC_SQRT normalization vector...
  [2023-09-29 17:12:19.123] [info]: [1000] DONE! Processed 26658348 pixels across 8 chromosomes in 10.07s
  ...
  [2023-09-29 17:12:37.412] [info]: DONE! Processed 10 resolution(s) in 28.43s!
  [2023-09-29 17:12:37.412] [info]: data/4DNFIZ1ZVXC8.hic9 size: 133.68 MB
  [2023-09-29 17:12:37.412] [info]: 4DNFIZ1ZVXC8.mcool size: 100.00 MB


It is also possible to convert only a subset of available resolutions by specifying resolutions to be converted with the ``--resolutions`` option.

When specifying a single resolution, the resulting file will be in .cool format.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.hic9 4DNFIZ1ZVXC8.1000.cool --resolutions 1000

  [2023-09-29 17:42:47.917] [info]: Running hictk v0.0.2-f83f93e
  [2023-09-29 17:42:47.917] [info]: Converting data/4DNFIZ1ZVXC8.hic9 to 4DNFIZ1ZVXC8.cool (hic -> cool)...
  [2023-09-29 17:42:47.982] [info]: [1000] begin processing 1000bp matrix...
  [2023-09-29 17:42:49.982] [info]: [1000] processing chr2R:11267000-11268000 at 5005005 pixels/s (cache hit rate 93.05%)...
  [2023-09-29 17:42:52.339] [info]: [1000] processing chr3R:5672000-5673000 at 4242681 pixels/s (cache hit rate 92.66%)...
  [2023-09-29 17:42:54.071] [info]: [1000] processing SCALE normalization vector...
  [2023-09-29 17:42:54.109] [info]: [1000] processing VC normalization vector...
  [2023-09-29 17:42:54.150] [info]: [1000] processing VC_SQRT normalization vector...
  [2023-09-29 17:42:54.931] [info]: [1000] DONE! Processed 26658348 pixels across 8 chromosomes in 6.95s
  [2023-09-29 17:42:54.931] [info]: DONE! Processed 1 resolution(s) in 7.01s!
  [2023-09-29 17:42:54.931] [info]: data/4DNFIZ1ZVXC8.hic9 size: 133.68 MB
  [2023-09-29 17:42:54.931] [info]: 4DNFIZ1ZVXC8.cool size: 36.74 MB



Converting from .[m]cool to .hic
--------------------------------

``hictk convert`` can also be used to convert .[m]cool files to .hic format.

The conversion steps are similar to those carried out to convert .hic to .[m]cool

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.mcool 4DNFIZ1ZVXC8.hic

  [2024-01-23 17:19:34.045] [info]: Running hictk v0.0.6-570037c-dirty
  [2024-01-23 17:19:34.045] [info]: Converting 4DNFIZ1ZVXC8.mcool to 4DNFIZ1ZVXC8.hic (mcool -> hic)...
  [2024-01-23 17:19:37.808] [info]: ingesting pixels at 2700513 pixels/s...
  [2024-01-23 17:19:41.916] [info]: ingesting pixels at 2434275 pixels/s...
  [2024-01-23 17:19:48.685] [info]: ingesting pixels at 2500000 pixels/s...
  [2024-01-23 17:19:52.753] [info]: ingesting pixels at 2458815 pixels/s...
  [2024-01-23 17:19:59.034] [info]: ingesting pixels at 2805049 pixels/s...
  [2024-01-23 17:20:07.190] [info]: writing header at offset 0
  [2024-01-23 17:20:07.190] [info]: begin writing interaction blocks to file "4DNFIZ1ZVXC8.hic"...
  [2024-01-23 17:20:07.190] [info]: [1000 bp] writing pixels for chr2L:chr2L matrix at offset 248...
  [2024-01-23 17:20:07.595] [info]: [1000 bp] written 2676654 pixels for chr2L:chr2L matrix
  [2024-01-23 17:20:07.651] [info]: [5000 bp] writing pixels for chr2L:chr2L matrix at offset 4303035...
  [2024-01-23 17:20:08.257] [info]: [5000 bp] written 2676654 pixels for chr2L:chr2L matrix
  [2024-01-23 17:20:08.366] [info]: [10000 bp] writing pixels for chr2L:chr2L matrix at offset 9144982...
  [2024-01-23 17:20:08.821] [info]: [10000 bp] written 1433133 pixels for chr2L:chr2L matrix
  ...
  [2024-01-23 17:21:30.092] [info]: [5000 bp] computing expected vector density
  [2024-01-23 17:21:30.240] [info]: [5000 bp] computing expected vector density
  [2024-01-23 17:21:30.297] [info]: [1000 bp] computing expected vector density
  [2024-01-23 17:21:30.784] [info]: [5000 bp] computing expected vector density
  [2024-01-23 17:21:30.784] [info]: writing 50 normalized expected value vectors at offset 142822186...
  [2024-01-23 17:21:30.785] [info]: writing 400 normalization vectors at offset 143709792...
  [2024-01-23 17:21:30.839] [info]: DONE! Processed 10 resolution(s) in 116.79s!
  [2024-01-23 17:21:30.839] [info]: 4DNFIZ1ZVXC8.mcool size: 139.38 MB
  [2024-01-23 17:21:30.839] [info]: 4DNFIZ1ZVXC8.hic size: 147.52 MB

**Tips:**

* When converting large .[m]cool files to .hic, ``hictk`` may need to create large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space.
* When converting .[m]cool files to .hic certain conversion steps can be performed in parallel. To improve performance, please make sure to increase the number of processing threads with option ``--thread``.
