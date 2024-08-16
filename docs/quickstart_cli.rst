..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Quickstart (CLI)
################

First, install hictk with one of the methods listed in the :doc:`Installation <./installation>` section.

Next, verify that hictk was installed correctly with:

.. code-block:: console

  user@dev:/tmp$ hictk --version
  hictk-v1.0.0

Command line interface
======================

hictk CLI support performing common operations on .hic and .cool files directly from the shell.

Verifying file integrity
------------------------

.. code-block:: sh

  hictk validate interactions.cool --validate-index

  hictk validate interactions.hic

For more detailed examples refer to :doc:`File validation <./file_validation>`

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

For more detailed examples refer to :doc:`Reading interactions <./reading_interactions>`

Other operations
----------------

* :doc:`Format conversion <./format_conversion>`
* :doc:`Creating .cool and .hic files <./creating_cool_and_hic_files>`
* :doc:`Converting single-resolution files to multi-resolution <./creating_multires_files>`
* :doc:`Balancing Hi-C matrices <./balancing_matrices>`


API
===

Refer to :doc:`Quickstart (API) <./quickstart_api>`.
