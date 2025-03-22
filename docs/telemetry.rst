..
   Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Telemetry
#########

Starting with version vx.x.x of hictk we introduced support for telemetry collection.

This only applies when hictk is invoked from the CLI (i.e. not to libhictk).

This page outlines what information we are collecting and why.
Furthermore, we provide instructions on how telemetry collection can be disabled at execution and compile time.

What information is being collected
-----------------------------------

``hictk`` is instrumented to collect general information about ``hictk`` itself and the system where it is being run.

We do not collect any sensitive information that could be used to identify our users, the machine or environment where ``hictk`` is being run, the datasets processed by ``hictk``, or the parameters used to run ``hictk``.

These is the data we are collecting:

* Information on how ``hictk`` was compiled (i.e. compiler name, version, and build type).
* Information on the system where ``hictk`` is being run (i.e. operating system and processor architecture).
* Information about ``hictk`` itself (i.e. version of ``hictk`` and of the library used for telemetry collection).
* How ``hictk`` is being invoked (i.e. the subcommand and input/output format(s) where applicable).
* Information about ``hictk`` execution (i.e. when it was launched, how long the command took to finish, and whether the command terminated with an error).
* For the ``hictk dump`` subcommand, we are also collecting the name of the table that is being dumped (e.g. pixels or chroms).

The following table shows an example of the telemetry collected when running ``hictk dump``:

.. csv-table:: Telemetry information collected when running ``hictk dump``
  :file: ./assets/telemetry_table.tsv
  :header-rows: 1
  :delim: tab

Why are we collecting this information?
---------------------------------------

There are two main motivations behind our decision to start collecting telemetry data:

#. To get an idea of how big our user base is - this will help us, among other things, to secure funding to maintain ``hictk`` in the future.
#. To better understand which of the functionalities offered by ``hictk`` are most used by our users - we intend to use this information to help us decide which features we should focus our development efforts on.

How is telemetry information processed and stored
-------------------------------------------------

Telemetry is sent to an OpenTelemetry collector running on a virtual server hosted on the Norwegian Research and Education Cloud (`NREC <https://www.nrec.no/>`_).

The virtual server and collector are managed by us, and traffic between ``hictk`` and the collector is encrypted.

The collector processes incoming data continuously and forwards it to a dashboard for data analytics and a backup solution (both services are hosted in Europe).
Communication between the collector, dashboard, and backup site is also encrypted.
Data stored by the dashboard and backup site is encrypted at rest.

The analytics dashboard keeps telemetry data for up to 60 days, while the backup site is currently set up to store telemetry data indefinitely (although this may change in the future).

How to disable telemetry collection
-----------------------------------

We provide two mechanisms to disable telemetry.

#. Disabling telemetry at runtime: simply define the ``HICTK_NO_TELEMETRY`` environment variable before launching ``hictk`` (e.g. ``HICTK_NO_TELEMETRY=1 hictk dump matrix.cool``)
#. Disabling telemetry at compile time: this only applies if you are building hictk from source as outlined in :doc:`installation_src`.

   To completely disable telemetry support at compile time pass ``-DHICTK_ENABLE_TELEMETRY=OFF`` when configuring the project with cmake.

   When ``HICTK_ENABLE_TELEMETRY`` is set to ``OFF``, classes and functions used to collect information using OpenTelemetry are replaced with alternative implementations that do nothing.
   Furthermore, the OpenTelemetry library is not linked to the ``hictk`` binary, meaning that no code involved in the collection of telemetry information is contained in or loaded by the ``hictk`` binary.

Where can I find the code used for telemetry collection?
--------------------------------------------------------

All code concerning telemetry collection is defined in the header file `hictk/tools/telemetry.hpp <https://github.com/paulsengroup/hictk/blob/main/src/hictk/include/hictk/tools/telemetry.hpp>`_.

The link flags and pre-processor macros toggling telemetry support at compile time are defined in `src/hictk/CMakeLists.txt <https://github.com/paulsengroup/hictk/blob/main/src/hictk/CMakeLists.txt>`_.
