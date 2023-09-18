..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Quickstart
##########

First, install hictk with one of the methods listed in the :doc:`Installation <./installation>` section.

Next, verify that hictk was installed correctly with:

.. code-block:: console

  user@dev:/tmp$ hictk --version
  hictk-v0.0.2

Command line interface
======================

hictk CLI support performing common operations on .hic and .cool files directly from the shell.

Verifying file integrity
------------------------

.. code-block:: sh

  hictk validate interactions.cool --validate-index

  hictk validate interactions.hic


Converting between file formats
-------------------------------

hictk supports converting matrices in .hic format to .cool/.mcool and viceversa using ``hictk convert``

Converting hic to cooler
^^^^^^^^^^^^^^^^^^^^^^^^

hictk has native support for converting .hic files to .cool and .mcool format:

.. code-block:: sh

  # Create a .mcool file using all resolutions available in interactions.hic
  hictk convert interactions.hic interactions.mcool

  # Create a .cool file at 10kb resolution
  hictk convert interactions.hic interactions.cool --resolutions 10000

  # Create a .mcool file using a subset of the resolutions available in interactions.hic
  hictk convert interactions.hic interactions.mcool --resolutions 10000 20000 50000

Converting cooler to hic
^^^^^^^^^^^^^^^^^^^^^^^^

hictk support converting .cool and .mcool files to .hic by wrapping `JuicerTools <https://github.com/aidenlab/juicertools>`_ or `HiCTools <https://github.com/aidenlab/HiCTools>`_.

First, download the .jar file for ``JuicerTools`` or ``HiCTools`` from e.g. `here <https://github.com/aidenlab/HiCTools/releases/latest>`_.

Next, test that ``JuicerTools`` or ``HiCTools`` can execute on your machine (this will require Java to be installed. Refer to ``JuicerTools`` documentation for more details.

.. code-block:: console

  user@dev:/tmp$ java -jar hic_tools.3.30.00.jar -help
  Juicer Tools Version 3.30.00
  Usage:
  	pre [options] <infile> <outfile> <genomeID>
  	addNorm <input_HiC_file> [input_vector_file]
  	pearsons [-p] <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
  	eigenvector -p <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr> <BP/FRAG> <binsize> [outfile]
  	sum [options] <outfile.hic> <infile1.hic> <infile2.hic> ... <infileN.hic
  	-h, --help print help
  	-v, --verbose verbose mode
  	-V, --version print version
  Type juicer_tools <commandName> for more detailed usage instructions


.. code-block:: sh

  # Create a .hic file using interactions.cool as base resolution
  hictk convert interactions.cool interactions.hic --juicer-tools-jar hic_tools.3.30.00.jar

  # Create a .hic file with the resolutions found in interactions.mcool
  hictk convert interactions.mcool interactions.cool --juicer-tools-jar hic_tools.3.30.00.jar

If you are getting an error like ``Exception in thread "main" java.lang.OutOfMemoryError: Java heap space`` this signals that ``JuicerTools`` is running out of memory.
In order to convert the matrix you should increase ``JuicerTools`` memory with e.g. ``--juicer-tools-memory=100G``.


Reading interactions
--------------------

hictk supports reading interactions from .hic and .cool files through the ``hictk dump`` command:

.. code-block:: console

  user@dev:/tmp$ hictk dump interactions.cool
  0	0	1745
  0	1	2844
  0	2	409
  ...

  user@dev:/tmp$ hictk dump interactions.cool --join
  chr2L	0	10000	chr2L	0	10000	1745
  chr2L	0	10000	chr2L	10000	20000	2844
  chr2L	0	10000	chr2L	20000	30000	409
  ...

  user@dev:/tmp$ hictk dump interactions.mcool::/resolutions/10000 --join
  chr2L	0	10000	chr2L	0	10000	1745
  chr2L	0	10000	chr2L	10000	20000	2844
  chr2L	0	10000	chr2L	20000	30000	409
  ...

  user@dev:/tmp$ hictk dump interactions.hic --join --resolution 10000 --matrix-type expected
  chr2L	0	10000	chr2L	0	10000	2351.23291015625
  chr2L	0	10000	chr2L	10000	20000	1447.001708984375
  chr2L	0	10000	chr2L	20000	30000	613.9473876953125
  ...

  user@dev:/tmp$ hictk dump interactions.hic --join --resolution 10000 --normalization VC
  chr2L	0	10000	chr2L	0	10000	3575.918701171875
  chr2L	0	10000	chr2L	10000	20000	2654.79052734375
  chr2L	0	10000	chr2L	20000	30000	387.9197082519531
  ...

  user@dev:/tmp$ hictk dump interactions.hic --join --resolution 10000 --range chr3L:20,000,000-25,000,000
  chr3L	20000000	20010000	chr3L	20000000	20010000	5400
  chr3L	20000000	20010000	chr3L	20010000	20020000	3766
  chr3L	20000000	20010000	chr3L	20020000	20030000	2015

  user@dev:/tmp$ hictk dump interactions.hic --join --resolution 10000 --range chr3L:20,000,000-25,000,000 --range2 chrX
  chr3L	20000000	20010000	chrX	50000	60000	2
  chr3L	20000000	20010000	chrX	140000	150000	1
  chr3L	20000000	20010000	chrX	150000	160000	1
  ...

Creating .cool files
--------------------

hictk supports creating .cool files from text files in the following formats:

* `pairs (4DN-DCIC) <https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md#example-pairs-file>`_
* `validPairs (nf-core/hic) <https://nf-co.re/hic/2.1.0/docs/output/#valid-pairs-detection-with-hic-pro>`_
* `bedGraph2 <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_
* `COO <https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays>`_

.. code-block:: sh

  # Create a 10kbp .cool file using hg38 as reference
  hictk load --format 4dn --assembly hg38 hg38.chrom.sizes 10000 out.cool < interactions.txt

  # Same as above but using gzip-compressed interactions
  zcat interactions.txt.gz | hictk load --format 4dn --assembly hg38 hg38.chrom.sizes 10000 out.cool

  # Using interactions in bedgraph2 format
  hictk load --format bg2 --assembly hg38 hg38.chrom.sizes 10000 out.cool < interactions.txt

Merging multiple Cooler files
-----------------------------

Multiple .cool files using the same reference genome and resolution can be merged using ``hictk merge``:

.. code-block:: sh

  hictk merge interactions1.cool interactions2.cool -o merged.cool

Converting .cool to .mcool
--------------------------

Interactions from a single-resolution Cooler file (.cool) can be used to generate a multi-resolution Cooler (.mcool) by iterative coaresning using ``hictk zoomify``

.. code-block:: sh

  hictk zoomify interactions.cool interactions.mcool --resolutions 1000 5000 10000 ...

  # Coarsen a single resolution
  hictk zoomify interactions.cool interactions.10000.cool --no-copy-base-resolution --resolutions 10000
