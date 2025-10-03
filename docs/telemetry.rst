..
   Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Telemetry
#########

Starting with version v2.1.0 of hictk, we introduced support for telemetry collection.

This only applies when hictk is invoked from the CLI (i.e., not to libhictk).

This page outlines what information we are collecting and why.
Furthermore, we provide instructions on how telemetry collection can be disabled at execution and compile time.

What information is being collected
-----------------------------------

``hictk`` is instrumented to collect general information about ``hictk`` itself and the system where it is being run.

We do not collect any sensitive information that could be used to identify our users, the machine or environment where ``hictk`` is being run, the datasets processed by ``hictk``, or the parameters used to run ``hictk``.

This is the data we are collecting:

* Information on how ``hictk`` was compiled (i.e., compiler name, version, and build type).
* Information on the system where ``hictk`` is being run (i.e., operating system and processor architecture).
* Information about ``hictk`` itself (i.e., version of ``hictk`` and its third-party version of dependencies).
* The continent, country, and region names, as well as the time zone where ``hictk`` was launched.
  This information is inferred from the IP address used to submit the telemetry (the IP address itself is not part of the telemetry data we collect and it never stored by our servers).
* How ``hictk`` is being invoked (i.e., the subcommand, the hash of the command line arguments used to invoke ``hictk``, and the input/output format(s) where applicable).
* Information about ``hictk`` execution (i.e., when it was launched, how long the command took to finish, and whether the command terminated with an error).
* For the ``hictk dump`` subcommand, we are also collecting the name of the table that is being dumped (e.g., pixels or chroms).

This is an example of the telemetry collected when running ``hictk dump``:

.. code-block:: text

  name          : subcommand.dump
  trace_id      : a51ce70f8aff91281eb70332c5eb775b
  span_id       : 869ec9f57d2e170e
  tracestate    :
  parent_span_id: 0000000000000000
  start         : 1758032386754347017
  duration      : 151590958
  description   :
  span kind     : Internal
  status        : Ok
  attributes    :
	param.table: pixels
	meta.input-format: mcool
	meta.argv-sha3-256: 6840f26c9293a323369ea6d571a48b8a49934e76e6e1a645f224caf14663c
	schema: 1
  events        :
  links         :
  resources     :
	build.compiler.name: Clang
	build.compiler.version: 20.1.8
	build.dependencies.boost.version: 1.88.0
	build.dependencies.bshoshany-thread-pool.version: 5.0.0
	build.dependencies.cli11.version: 2.5.0
	build.dependencies.concurrentqueue.version: 1.0.4
	build.dependencies.fast_float.version: 8.0.2
	build.dependencies.fmt.version: 11.2.0
	build.dependencies.hdf5.version: 1.14.6
	build.dependencies.highfive.version: 2.10.0
	build.dependencies.libarchive.version: 3.8.1
	build.dependencies.libdeflate.version: 1.23
	build.dependencies.nlohmann_json.version: 3.12.0
	build.dependencies.opentelemetry-cpp.version: 1.21.0
	build.dependencies.parallel-hashmap.version: 2.0.0
	build.dependencies.readerwriterqueue.version: 1.0.6
	build.dependencies.span-lite.version: 0.11.0
	build.dependencies.tomlplusplus.version: 3.4.0
	build.dependencies.zstd.version: 1.5.7
	build.type: Release
  geo.continent_name: Europe
  geo.country_name: Norway
  geo.region_name: Oslo County
  geo.timezone: Europe/Oslo
	host.arch: x86_64
	os.type: Linux
	os.version: 6.16.3-200.fc42.x86_64
	service.name: hictk
	service.version: 2.1.5
	telemetry.sdk.language: cpp
	telemetry.sdk.name: opentelemetry
	telemetry.sdk.version: 1.21.0
  instr-lib     : hictk

Why are we collecting this information?
---------------------------------------

There are two main motivations behind our decision to start collecting telemetry data:

#. To get an idea of how big our user base is: this will help us, among other things, to secure funding to maintain ``hictk`` in the future.
#. To better understand which of the functionalities offered by ``hictk`` are most used by our users: we intend to use this information to help us decide which features we should focus our development efforts on.

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

#. Disabling telemetry at runtime: simply define the ``HICTK_NO_TELEMETRY`` environment variable before launching ``hictk`` (e.g., ``HICTK_NO_TELEMETRY=1 hictk dump matrix.cool``)
#. Disabling telemetry at compile time: this only applies if you are building hictk from source as outlined in :doc:`installation_src`.

   To completely disable telemetry support at compile time pass ``-DHICTK_ENABLE_TELEMETRY=OFF`` when configuring the project with CMake.

   When ``HICTK_ENABLE_TELEMETRY`` is set to ``OFF``, classes and functions used to collect information using OpenTelemetry are replaced with alternative implementations that do nothing.
   Furthermore, the OpenTelemetry library is not linked to the ``hictk`` binary, meaning that no code involved in the collection of telemetry information is contained in or loaded by the ``hictk`` binary.

Where can I find the code used for telemetry collection?
--------------------------------------------------------

All code concerning telemetry collection is defined in the library under `src/hictk/telemetry <https://github.com/paulsengroup/hictk/tree/main/src/hictk/telemetry>`_.

The link flags and pre-processor macros toggling telemetry support at compile time are defined in files `src/hictk/CMakeLists.txt <https://github.com/paulsengroup/hictk/blob/main/src/hictk/CMakeLists.txt>`_ and `src/hictk/telemetry/CMakeLists.txt <https://github.com/paulsengroup/hictk/blob/main/src/hictk/telemetry/CMakeLists.txt>`_.
