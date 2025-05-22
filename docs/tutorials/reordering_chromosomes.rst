..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Reorder chromosomes
###################

TLDR
----

.. code-block:: console

    # Important! --bin-size should be the same resolution as matrix.cool
    user@dev:/tmp hictk load <(hictk dump --join matrix.cool) \
                             output.cool \
                             --chrom-sizes=<(hictk dump --table=chroms matrix.cool | sort -k2,2nr) \
                             --format=bg2 \
                             --bin-size=1kbp \
                             --transpose-lower-triangular-pixels


Why is this needed?
-------------------

Sometimes we want to compare files using the same reference genome assembly, but with different chromosome orders (e.g. in one file chromosomes are sorted by size while in the other they are sorted by name).
This can be a problem especially when trying to visually compare such files.
This tutorial shows how to convert a .cool file with chromosomes sorted by name to a .cool file with chromosomes sorted by size.
The same procedure can be applied to .hic files.


Walkthrough
-----------

For this tutorial, we will use file ``4DNFIOTPSS3L.hic`` as an example, which can be downloaded from `here <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic>`__.

First, we extract the list of chromosomes from the input file:

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIOTPSS3L.hic --table=chroms | tee chrom.sizes

    2L	23513712
    2R	25286936
    3L	28110227
    3R	32079331
    4	1348131
    X	23542271
    Y	3667352

Second, we re-order chromosomes:

.. code-block:: console

    user@dev:/tmp sort -k2,2nr chrom.sizes | tee chrom.sizes.sorted

    3R	32079331
    3L	28110227
    2R	25286936
    X	23542271
    2L	23513712
    Y	3667352
    4	1348131

Next, we dump pixels in bedGraph2 format (see below for how to make this step more efficient):

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIOTPSS3L.hic --join --resolution 1kbp > pixels.bg2

    user@dev:/tmp head pixels.bg2

    2L	5000	6000	2L	5000	6000	41
    2L	5000	6000	2L	6000	7000	126
    2L	5000	6000	2L	7000	8000	60
    2L	5000	6000	2L	8000	9000	77
    2L	5000	6000	2L	9000	10000	97
    2L	5000	6000	2L	10000	11000	3
    2L	5000	6000	2L	11000	12000	1
    2L	5000	6000	2L	12000	13000	66
    2L	5000	6000	2L	13000	14000	116
    2L	5000	6000	2L	14000	15000	64

Finally, we load pixels into a new .hic file

.. code-block:: console

    user@dev:/tmp hictk load pixels.bg2 \
                             output.hic \
                             --chrom-sizes=chrom.sizes.sorted \
                             --transpose-lower-triangular-pixels \
                             --format=bg2 \
                             --bin-size=1kbp

    [2024-09-27 19:00:40.344] [info]: Running hictk v1.0.0-fbdcb591
    [2024-09-27 19:00:40.353] [info]: begin loading pixels into a .hic file...
    [2024-09-27 19:00:42.504] [info]: preprocessing chunk #1 at 4847310 pixels/s...
    [2024-09-27 19:00:45.244] [info]: preprocessing chunk #2 at 3649635 pixels/s...
    [2024-09-27 19:00:48.180] [info]: preprocessing chunk #3 at 3407155 pixels/s...
    [2024-09-27 19:00:50.616] [info]: preprocessing chunk #4 at 4105090 pixels/s...
    [2024-09-27 19:00:53.251] [info]: preprocessing chunk #5 at 3203434 pixels/s...
    [2024-09-27 19:00:54.358] [info]: writing header at offset 0
    [2024-09-27 19:00:54.358] [info]: begin writing interaction blocks to file "output.hic"...
    [2024-09-27 19:00:54.358] [info]: [1000 bp] writing pixels for 3R:3R matrix at offset 171...
    [2024-09-27 19:01:01.039] [info]: [1000 bp] written 9571521 pixels for 3R:3R matrix
    ...
    [2024-09-27 19:01:26.831] [info]: [1000 bp] initializing expected value vector
    [2024-09-27 19:01:32.649] [info]: [1000 bp] computing expected vector density
    [2024-09-27 19:01:32.649] [info]: writing 1 expected value vectors at offset 93720080...
    [2024-09-27 19:01:32.649] [info]: writing 0 normalized expected value vectors at offset 93848475...
    [2024-09-27 19:01:32.682] [info]: ingested 114355295 interactions (48437845 nnz) in 52.337885908s!

Lastly, we check that chromosomes are properly sorted:

.. code-block:: console

    user@dev:/tmp hictk dump output.hic --table=chroms

    3R	32079331
    3L	28110227
    2R	25286936
    X	23542271
    2L	23513712
    Y	3667352
    4	1348131

Tips and tricks
---------------

There is one potential problem with the above solution, and that is the size of file ``pixels.bg2``
Luckily, we can completely avoid generating this file by using output redirection and process substitutions:

.. code-block:: console

    user@dev:/tmp hictk load <(hictk dump 4DNFIOTPSS3L.hic --join --resolution 1kbp) \
                             output.hic \
                             --chrom-sizes=chrom.sizes.sorted \
                             --transpose-lower-triangular-pixels \
                             --format=bg2 \
                             --bin-size=1kbp


Note that hictk still needs to generate some temporary file to load interactions into a new .cool or .hic file.
When processing large files, it is a good idea to specify custom folder where to create temporary files through the ``--tmpdir`` flag:

.. code-block:: console

    user@dev:/tmp hictk load <(hictk dump 4DNFIOTPSS3L.hic --join --resolution 1kbp) \
                             output.hic \
                             --chrom-sizes=chrom.sizes.sorted \
                             --transpose-lower-triangular-pixels \
                             --format=bg2 \
                             --bin-size=1kbp \
                             --tmpdir=/var/tmp/

Another option you may want to consider when working with .hic files is the ``--threads`` option, which can significantly reduce the time required to load interactions into .hic files.
