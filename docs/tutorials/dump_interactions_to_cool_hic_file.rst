..
   Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Dump interactions to .cool or .hic file
#######################################

TLDR
----

.. code-block:: console

    # Important! --bin-size should be the same resolution as matrix.cool
    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             <(hictk dump --table=chroms matrix.cool)
                             output.cool \
                             < <(hictk dump --join
                                            --range=chr2L:0-10,000,000
                                            --range2=chr3R:0-10,000,000
                                            matrix.cool)

Why is this needed?
-------------------

``hictk dump`` can read interactions from .cool, .mcool, and .hic files and write them in text format to stdout.
Additionally, ``hictk dump`` supports fetching interactions overlapping a pair of regions of interest through the ``--range`` and ``--range2`` CLI options.
However, instead of writing interactions to stdout, we may want to write them to a new .cool or .hic file.
This tutorial shows how this can be accomplished using ``hictk dump`` and ``hictk load``.


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

Second, we dump pixels in bedGraph2 format (see below for how to make this step more efficient):

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIZ1ZVXC8.mcool \
                             --join \
                             --resolution=1000 \
                             --range=chr2L:5,000,000-10,000,000 \
                             --range2=chr3R:7,500,000-10,000,000 > pixels.bg2

    user@dev:/tmp head pixels.bg2

    chr2L	5000000	5001000	chr3R	7506000	7507000	1
    chr2L	5000000	5001000	chr3R	7624000	7625000	1
    chr2L	5000000	5001000	chr3R	7943000	7944000	1
    chr2L	5000000	5001000	chr3R	8014000	8015000	1
    chr2L	5000000	5001000	chr3R	8130000	8131000	1
    chr2L	5000000	5001000	chr3R	8245000	8246000	1
    chr2L	5000000	5001000	chr3R	8855000	8856000	1
    chr2L	5000000	5001000	chr3R	9032000	9033000	1
    chr2L	5000000	5001000	chr3R	9171000	9172000	1
    chr2L	5000000	5001000	chr3R	9380000	9381000	1


Finally, we load pixels into a new .cool file

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             chrom.sizes \
                             output.cool < pixels.bg2

    [2024-03-21 13:22:57.542] [info]: Running hictk v0.0.10-1c2bafd
    [2024-03-21 13:22:57.542] [info]: begin loading unsorted pixels into a .cool file...
    [2024-03-21 13:22:57.613] [info]: writing chunk #1 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 13:22:57.630] [info]: done writing chunk #1 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 13:22:57.630] [info]: writing chunk #2 to intermediate file "/tmp/output.cool.tmp/output.cool.tmp"...
    [2024-03-21 13:22:57.634] [info]: done writing chunk #2 to tmp file "/tmp/output.cool.tmp/output.cool.tmp".
    [2024-03-21 13:22:57.634] [info]: merging 2 chunks into "output.cool"...
    [2024-03-21 13:22:57.676] [info]: ingested 26214 interactions (25085 nnz) in 0.133955616s!

Removing empty chromosomes from the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This can be easily achieved by grepping ``chr2L`` and ``chr3R`` when generating the ``chrom.sizes`` file.

.. code-block:: console

    user@dev:/tmp hictk dump 4DNFIZ1ZVXC8.mcool --table=chroms |
                  grep -e 'chr2L' -e 'chr3R' |
                  tee chrom.sizes

    chr2L	23513712
    chr3R	32079331


Tips and tricks
---------------

There is one potential problem with the above solution, and that is the size of file ``pixels.bg2``
Luckily, we can completely avoid generating this file by using output redirection and process substitutions:

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             chrom.sizes \
                             output.cool \
                             < <(hictk dump 4DNFIZ1ZVXC8.mcool \
                                            --join \
                                            --resolution=1000 \
                                            --range=chr2L:0-10,000,000 \
                                            --range2=chr3R:0-10,000,000)

Note that hictk still needs to generate some temporary file to load interactions into a new .cool or .hic file.
When processing large files, it is a good idea to specify custom folder where to create temporary files through the ``--tmpdir`` flag:

.. code-block:: console

    user@dev:/tmp hictk load --format=bg2 \
                             --bin-size=1000 \
                             --tmpdir=/var/tmp/ \
                             chrom.sizes.sorted \
                             output.cool \
                             < <(hictk dump 4DNFIZ1ZVXC8.mcool \
                                            --join \
                                            --resolution=1000 \
                                            --range=chr2L:0-10,000,000 \
                                            --range2=chr3R:0-10,000,000)

Another option you may want to consider when working with .hic files, is the ``--threads`` option, which can significantly reduce the time required to load interactions into .hic files.
