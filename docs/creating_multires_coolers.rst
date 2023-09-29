..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Converting .cool to .mcool
##########################

Interactions from a single-resolution Cooler file (.cool) can be used to generate a multi-resolution Cooler (.mcool) by iterative coaresning using ``hictk zoomify``

.. code-block:: console

  user@dev:/tmp$ hictk zoomify data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 out.mcool

  [2023-09-29 19:28:39.926] [info]: Running hictk v0.0.2
  [2023-09-29 19:28:39.929] [info]: coarsening cooler at data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 13 times (1000 -> 1000 -> 2000 -> 5000 -> 10000 -> 20000 -> 50000 -> 100000 -> 200000 -> 500000 -> 1000000 -> 2000000 -> 5000000 -> 10000000)
  [2023-09-29 19:28:39.929] [info]: copying 1000 resolution from data/4DNFIZ1ZVXC8.mcool::/resolutions/1000
  [2023-09-29 19:28:40.119] [info]: generating 2000 resolution from 1000 (2x)
  [2023-09-29 19:28:40.343] [info]: [1000 -> 2000] processing chr2L:1996000-1998000 at 4484305 pixels/s...
  [2023-09-29 19:28:40.663] [info]: [1000 -> 2000] processing chr2L:4932000-4934000 at 3125000 pixels/s...
  [2023-09-29 19:28:40.973] [info]: [1000 -> 2000] processing chr2L:7986000-7988000 at 3236246 pixels/s...
  ...
  [2023-09-29 19:29:12.513] [info]: generating 10000000 resolution from 5000000 (2x)
  [2023-09-29 19:29:12.519] [info]: DONE! Processed 13 resolution(s) in 32.59s!

  # Coarsen a single resolution
  user@dev:/tmp$ hictk zoomify data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 out.cool --resolutions 50000

  [2023-09-29 19:30:52.476] [info]: Running hictk v0.0.2
  [2023-09-29 19:30:52.482] [info]: coarsening cooler at data/4DNFIZ1ZVXC8.mcool::/resolutions/1000 2 times (1000 -> 1000 -> 50000)
  [2023-09-29 19:30:52.482] [info]: copying 1000 resolution from data/4DNFIZ1ZVXC8.mcool::/resolutions/1000
  [2023-09-29 19:30:52.668] [info]: generating 50000 resolution from 1000 (50x)
  [2023-09-29 19:30:53.789] [info]: [1000 -> 50000] processing chr2L:23000000-23050000 at 896057 pixels/s...
  [2023-09-29 19:30:55.005] [info]: [1000 -> 50000] processing chr3L:4600000-4650000 at 822368 pixels/s...
  [2023-09-29 19:30:56.440] [info]: [1000 -> 50000] processing chr3R:32050000-32079331 at 696864 pixels/s...
  [2023-09-29 19:30:56.863] [info]: DONE! Processed 2 resolution(s) in 4.39s!
