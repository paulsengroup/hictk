..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Balancing Hi-C matrices
#######################

``hictk`` supports balancing .hic, .cool and .mcool files using ICE (iterative correction and eigenvector decomposition), SCALE and VC:

.. code-block:: console

  user@dev:/tmp$ hictk balance --help
  Balance Hi-C matrices using ICE, SCALE, or VC.
  Usage: hictk balance [OPTIONS] [SUBCOMMAND]

  Options:
    -h,--help                   Print this help message and exit

  Subcommands:
    ice                         Balance Hi-C matrices using ICE.
    scale                       Balance Hi-C matrices using SCALE.
    vc                          Balance Hi-C matrices using VC.

The following is an example showing how to balance a .cool file using ICE.

.. code-block:: console

  user@dev:/tmp$ hictk balance ice 4DNFIZ1ZVXC8.mcool::/resolutions/1000

  [2024-09-26 16:02:19.731] [info]: Running hictk v1.0.0-fbdcb591
  [2024-09-26 16:02:19.731] [info]: balancing using ICE (GW_ICE)
  [2024-09-26 16:02:19.734] [info]: Writing interactions to temporary file /tmp/hictk-tmp-XXXX1ZC9FF/4DNFIZ1ZVXC8.mcool.tmp...
  [2024-09-26 16:02:22.480] [info]: Initializing bias vector...
  [2024-09-26 16:02:22.482] [info]: Masking rows with fewer than 10 nnz entries...
  [2024-09-26 16:02:23.392] [info]: Masking rows using mad_max=5...
  [2024-09-26 16:02:23.860] [info]: Iteration 1: 36452362.243888594
  [2024-09-26 16:02:24.327] [info]: Iteration 2: 21649057.88060747
  [2024-09-26 16:02:24.792] [info]: Iteration 3: 7890065.688497526
  ...
  [2024-09-26 16:03:12.285] [info]: Iteration 107: 2.0533518142916073e-05
  [2024-09-26 16:03:12.752] [info]: Iteration 108: 1.601698258037195e-05
  [2024-09-26 16:03:13.216] [info]: Iteration 109: 1.2493901433163442e-05
  [2024-09-26 16:03:13.681] [info]: Iteration 110: 9.745791018854495e-06
  [2024-09-26 16:03:13.707] [info]: Writing weights to 4DNFIZ1ZVXC8.mcool::/resolutions/1000/bins/GW_ICE...
  [2024-09-26 16:03:13.708] [info]: Linking weights to 4DNFIZ1ZVXC8.mcool::/resolutions/1000/bins/weight...

When balancing files in .mcool or .hic formats, all resolutions are balanced.

By default balancing coefficients are stored in the input file under the name of "weight".

This can be changed by passing the desired name through the ``--name`` option.

``hictk`` supports three balancing methods:

* Using all (genome-wide) interactions (default)
* Using trans interactions only
* Using cis interactions only

Balancing method can be changed through the ``--mode`` option (e.g. ``--mode=gw`` or ``--mode=cis``).

When enough memory is available, ``hictk`` can be instructed to load all interactions into system memory by passing the ``--in-memory`` flag. This can dramatically speed up matrix balancing at the cost of potentially much higher memory usage (approximately 1 GB of RAM for every 40M interactions).

Another way to improve performance is to increase the number of threads available for computation using the ``--thread`` option.
It should be noted that when using a large number of threads (e.g. more than 16) without the ``--in-memory`` option, performance is likely limited by disk throughput. Thus, users are advised to use a large number of threads only when temporary data (``/tmp`` by default on most UNIX-like systems) is stored on a fast SSD.

When the ``--in-memory`` option is not used, ``hictk`` will create a temporary file under the default temporary folder. This file stores interactions using a layout and compression that are optimized for the access pattern used by ``hictk balance``. When balancing large matrices, this file can be quite large (sometimes tens of GBs). If this is the case, it may be appropriate to change the temporary folder using the ``--tmpdir`` option.

Finally, when balancing .hic files, only .hic v9 files and newer are supported.
