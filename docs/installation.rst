..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Installation
############

Conda (bioconda)
================

The hictk package for Linux and macOS is available on bioconda and can be installed as follows:

.. code-block:: console

  user@dev:/tmp$ conda create -n hictk -c conda-forge -c bioconda hictk

  user@dev:/tmp$ conda activate hictk

  (hictk) user@dev:/tmp$ whereis hictk
  hictk: /home/user/.miniconda3/envs/hictk/bin/hictk

  (hictk) user@dev:/tmp$ hictk --version
  hictk-v2.1.4-bioconda


Containers (Docker or Singularity/Apptainer)
============================================

First, ensure you have followed the instructions on how to install Docker or Singularity/Apptainer on your OS.

.. raw:: html

   <details>
   <summary><a>Installing Docker</a></summary>

The following instructions assume you have root/admin permissions.

* `Linux <https://docs.docker.com/desktop/install/linux-install/>`_
* `macOS <https://docs.docker.com/desktop/install/mac-install/>`_
* `Windows <https://docs.docker.com/desktop/install/windows-install/>`_

On some Linux distributions, simply installing Docker is not enough.
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
and `DockerHub <https://hub.docker.com/r/paulsengroup/hictk>`_.

Downloading and running the latest stable release can be done as follows:

.. code-block:: console

  # Using Docker, may require sudo
  user@dev:/tmp$ docker run ghcr.io/paulsengroup/hictk:2.1.4 --help

  # Using Singularity/Apptainer
  user@dev:/tmp$ singularity run ghcr.io/paulsengroup/hictk:2.1.4 --help

  Blazing fast tools to work with .hic and .cool files.
  Usage: hictk [OPTIONS] SUBCOMMAND
  Options:
    -h,--help                   Print this help message and exit
    -V,--version                Display program version information and exit
  Subcommands:
    balance                     Balance Hi-C files using ICE, SCALE, or VC.
    convert                     Convert Hi-C files between different formats.
    dump                        Read interactions and other kinds of data from .hic and Cooler files and write them to stdout.
    fix-mcool                   Fix corrupted .mcool files.
    load                        Build .cool and .hic files from interactions in various text formats.
    merge                       Merge multiple Cooler or .hic files into a single file.
    metadata                    Print file metadata to stdout.
    rename-chromosomes, rename-chroms
                                Rename chromosomes found in Cooler files.
    validate                    Validate .hic and Cooler files.
    zoomify                     Convert single-resolution Cooler and .hic files to multi-resolution by coarsening.

The above will print hictk's help message, and is equivalent to running :code:`hictk --help` from the command line (assuming hictk is available on your machine).

Installing from source
======================

Please refer to hictk's :doc:`build instructions <./installation_src>`.
