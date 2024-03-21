..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Reordering chromosomes
######################

TLDR
----

.. code-block:: console

    # Important! --bin-size should be the same resolution as matrix.cool
    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             <(hictk dump --table=chroms matrix.cool |
                               sort -k2,2nr) \
                             output.cool \
                             < <(hictk dump --join matrix.cool)


Why is this needed?
-------------------

Sometimes we want to compare files using the same reference genome assembly, but with different chromosome orders (e.g. in one file chromosomes are sorted by size while in the other they are sorted by name).
This can be a problem especially when trying to visually compare such files.
This tutorial shows how to convert a .cool file with chromosomes sorted by name to a .cool file with chromosomes sorted by size.
The same procedure can be applied to .hic files.


Walkthrough
-----------

For this tutorial, we will use file ``4DNFIZ1ZVXC8.mcool`` as an example, which can be downloaded from `here <https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1cf3518f-839a-42b9-b2c7-7f81ad5935c3/4DNFIZ1ZVXC8.mcool>`__.

First, we extract the list of chromosomes from the input file:

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIZ1ZVXC8.mcool --table=chroms | tee chrom.sizes

    chr2L	23513712
    chr2R	25286936
    chr3L	28110227
    chr3R	32079331
    chr4	1348131
    chrX	23542271
    chrY	3667352

Second, we re-order chromosomes:

.. code-block:: console

    user@dev:/tmp sort -k2,2nr chrom.sizes | tee chrom.sizes.sorted

    chr3R	32079331
    chr3L	28110227
    chr2R	25286936
    chrX	23542271
    chr2L	23513712
    chrY	3667352
    chr4	1348131

Next, we dump pixels in bedGraph2 format (see below for how to make this step more efficient):

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIZ1ZVXC8.mcool --join --resolution=1000 > pixels.bg2

    user@dev:/tmp head pixels.bg2

    chr2L	5000	6000	chr2L	5000	6000	127
    chr2L	5000	6000	chr2L	6000	7000	129
    chr2L	5000	6000	chr2L	7000	8000	60
    chr2L	5000	6000	chr2L	8000	9000	77
    chr2L	5000	6000	chr2L	9000	10000	97
    chr2L	5000	6000	chr2L	10000	11000	3
    chr2L	5000	6000	chr2L	11000	12000	1
    chr2L	5000	6000	chr2L	12000	13000	66
    chr2L	5000	6000	chr2L	13000	14000	116
    chr2L	5000	6000	chr2L	14000	15000	64

Finally, we load pixels into a new .cool file

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             chrom.sizes.sorted \
                             output.cool < pixels.bg2

    [2024-03-21 12:27:16.998] [info]: Running hictk v0.0.10-1c2bafd
    [2024-03-21 12:27:16.998] [info]: begin loading unsorted pixels into a .cool file...
    [2024-03-21 12:27:17.077] [info]: writing chunk #1 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:20.945] [info]: done writing chunk #1 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:20.945] [info]: writing chunk #2 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:24.890] [info]: done writing chunk #2 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:24.890] [info]: writing chunk #3 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:28.823] [info]: done writing chunk #3 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:28.823] [info]: writing chunk #4 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:32.668] [info]: done writing chunk #4 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:32.668] [info]: writing chunk #5 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:36.070] [info]: done writing chunk #5 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:36.070] [info]: writing chunk #6 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 12:27:36.079] [info]: done writing chunk #6 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 12:27:36.080] [info]: merging 6 chunks into "output.cool"...
    [2024-03-21 12:27:38.572] [info]: processing chr3R:20786000-20787000 chr3R:20808000-20809000 at 4091653 pixels/s...
    [2024-03-21 12:27:41.443] [info]: processing chr3L:7391000-7392000 chr3L:7417000-7418000 at 3484321 pixels/s...
    [2024-03-21 12:27:44.292] [info]: processing chr2R:9278000-9279000 chrX:5993000-5994000 at 3510004 pixels/s...
    [2024-03-21 12:27:47.062] [info]: processing chrX:14217000-14218000 chrX:17476000-17477000 at 3611412 pixels/s...
    [2024-03-21 12:27:49.901] [info]: ingested 119208613 interactions (48469783 nnz) in 32.902465965s!


Lastly, we check that chromosomes are properly sorted:

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIZ1ZVXC8.mcool --table=chroms

    chr3R	32079331
    chr3L	28110227
    chr2R	25286936
    chrX	23542271
    chr2L	23513712
    chrY	3667352
    chr4	1348131


Tips and tricks
---------------

There is one potential problem with the above solution, and that is the size of file ``pixels.bg2``
Luckily, we can completely avoid generating this file by using output redirection and process substitutions:

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             chrom.sizes.sorted \
                             output.cool \
                             < <(hictk dump 4DNFIZ1ZVXC8.mcool --join --resolution=1000)

Note that hictk still needs to generate some temporary file to load interactions into a new .cool or .hic file.
When processing large files, it is a good idea to specify custom folder where to create temporary files through the ``--tmpdir`` flag:

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             --tmpdir=/var/tmp/ \
                             chrom.sizes.sorted \
                             output.cool \
                             < <(hictk dump 4DNFIZ1ZVXC8.mcool --join --resolution=1000)

Another option you may want to consider when working with .hic files, is the ``--threads`` option, which can significantly reduce the time required to load interactions into .hic files.
