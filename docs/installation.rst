..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Installation
############

Conda (bioconda)
================

hictk package for Linux and MacOS is available on bioconda and can be installed as follows:

.. code-block:: console

  user@dev:/tmp$ conda create -n hictk -c conda-forge -c bioconda hictk

  user@dev:/tmp$ conda activate hictk

  (hictk) user@dev:/tmp$ whereis hictk
  hictk: /home/user/.miniconda3/envs/hictk/bin/hictk

  (hictk) user@dev:/tmp$ hictk --version
  hictk-v0.0.11-bioconda

If you are trying to install hictk on a Mac with an M chip, the above command may fail due to conda not being able to find a package for hictk.
You can workaround the above issue by prefixing conda commands with :code:`CONDA_SUBDIR=osx-64`.
Note that this will make hictk quite a bit slower, as the installed binary will be executed through Rosetta.
If performance is important, please consider :doc:`compiling hictk from source <./installation_src>` or using containers (see below).

Containers (Docker or Singularity/Apptainer)
============================================

First, make sure you follow the instructions on how to install Docker or Singularity/Apptainer on your OS.

.. raw:: html

   <details>
   <summary><a>Installing Docker</a></summary>

The following instructions assume you have root/admin permissions.

* `Linux <https://docs.docker.com/desktop/install/linux-install/#generic-installation-steps/>`_
* `MacOS <https://docs.docker.com/desktop/install/mac-install/>`_
* `Windows <https://docs.docker.com/desktop/install/windows-install/>`_

On some Linux distributions just installing Docker is not enough.
You also need to start (and optionally enable) the appropriate service(s).
This is usually done with one of the following:

.. code-block:: sh

  sudo systemctl start docker
  sudo systemctl start docker.service


Refer to `Docker <https://docs.docker.com/engine/install/>`_ or your OS/distribution documentation for more details.

.. raw:: html

   </details>

Pulling hictk Docker image
--------------------------

hictk Docker images are available on `GHCR.io <https://github.com/paulsengroup/hictk/pkgs/container/hictk>`_
and `DockerHub <https://hub.docker.com/repository/docker/paulsengroup/hictk>`_.

Downloading and running the latest stable release can be done as follows:

.. code-block:: console

  # Using Docker, may require sudo
  user@dev:/tmp$ docker run ghcr.io/paulsengroup/hictk:0.0.11 --help

  # Using Singularity/Apptainer
  user@dev:/tmp$ singularity run ghcr.io/paulsengroup/hictk:0.0.11 --help

  Blazing fast tools to work with .hic and .cool files.
  Usage: /usr/local/bin/hictk [OPTIONS] SUBCOMMAND

  Options:
    -h,--help                   Print this help message and exit
    -V,--version                Display program version information and exit

  Subcommands:
    convert                     Convert HiC matrices to a different format.
    dump                        Dump data from .hic and Cooler files to stdout.
    load                        Build .cool files from interactions in various text formats.
    merge                       Merge coolers.
    validate                    Validate .hic and Cooler files.
    zoomify                     Convert single-resolution Cooler file to multi-resolution by coarsening.

The above will print hictk's help message, and is equivalent to running :code:`hictk --help` on the command line (assuming hictk is available on your machine).

Installing from source
======================

Please refer to hictk's :doc:`build instructions <./installation_src>`.
