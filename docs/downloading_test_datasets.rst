..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Downloading test datasets
#########################

Test dataset for ``hictk`` are hosted on Zenodo: `doi.org/10.5281/zenodo.8121686 <https://doi.org/10.5281/zenodo.8121686>`_

After downloading the data, move to a folder with at least ~1 GB of free space and extract the test datasets:

.. code-block:: console
  :class: no-copybutton

  user@dev:/tmp$ mkdir data/
  user@dev:/tmp$ tar -xf hictk_test_data.tar.xz                     \
                     -C data --strip-components=3                   \
                     test/data/hic/4DNFIZ1ZVXC8.hic9                \
                     test/data/integration_tests/4DNFIZ1ZVXC8.mcool \
                     test/data/integration_tests/4DNFIKNWM36K.subset.pairs.xz

  user@dev:/tmp$ ls -lah data
  total 261M
  drwx------  2 dev dev   80 Sep 29 17:00 .
  drwxrwxrwt 26 dev dev  960 Sep 29 17:00 ..
  -rw-------  1 dev dev 128M Jun  8 19:42 4DNFIZ1ZVXC8.hic9
  -rw-------  1 dev dev 133M Jul  7 16:29 4DNFIZ1ZVXC8.mcool
