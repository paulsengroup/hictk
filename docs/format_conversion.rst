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

As ``libhictk`` is not yet capable of writing .hic files, ``hictk convert`` requires `JuicerTools <https://github.com/aidenlab/Juicebox/releases/latest>`_ or `HiCTools <https://github.com/aidenlab/HiCTools/releases/latest>`_ JARs and java to be available.

You should use HiCTools JAR unless you need to open the resulting .hic file with applications that do not support the latest .hic format specification.

.. code-block:: console

  user@dev:/tmp$ hictk convert data/4DNFIZ1ZVXC8.mcool 4DNFIZ1ZVXC8.hic --juicer-tools-jar hic_tools.3.30.00.jar

  [2023-09-29 17:44:10.001] [info]: Running hictk v0.0.2-f83f93e
  [2023-09-29 17:44:10.001] [info]: Converting data/4DNFIZ1ZVXC8.mcool to 4DNFIZ1ZVXC8.hic (mcool -> hic)...
  [2023-09-29 17:44:10.004] [info]: writing chromosomes to file /tmp/hictk-tmp-XXXXjjxVhi/reference.chrom.sizes...
  [2023-09-29 17:44:10.004] [info]: DONE! Wrote 8 chromosomes to file /tmp/hictk-tmp-XXXXjjxVhi/reference.chrom.sizes
  [2023-09-29 17:44:10.004] [info]: writing pixels to file /tmp/hictk-tmp-XXXXjjxVhi/pixels.tsv.gz...
  [2023-09-29 17:44:19.933] [info]: processing chr2R:19727000-19728000 chr2R:21162000-21163000 at 1007252 pixels/s...
  [2023-09-29 17:44:31.012] [info]: processing chr3R:6457000-6458000 chr3R:21482000-21483000 at 902609 pixels/s...
  [2023-09-29 17:44:37.397] [info]: wrote 26591454 pixels across 8 chromosomes to /tmp/hictk-tmp-XXXXjjxVhi/pixels.tsv.gz in 27.39s
  [2023-09-29 17:44:37.398] [info]: running juicer_tools pre...
  ...

**Tips:**

* Use JuicerTools instead of HiCTools if the output .hic file needs to be opened by applications that do not support the latest .hic format specification.
* When converting large .[m]cool files to .hic, ``hictk`` needs to create large temporary files. When this is the case, use option ``--tmpdir`` to set the temporary folder to a path with sufficient space
* When converting .[m]cool files to .hic, ``hictk`` tries to use ``pigz`` instead of plain ``gzip`` to compress temporary files. This can drammatically reduce conversion time. Please make sure ``pigz`` is installed and increase the number of processing thread with option ``--thread``.
