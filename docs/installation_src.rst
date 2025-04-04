..
   Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
   SPDX-License-Identifier: MIT

Installation (source)
#####################

Instructions assume hictk is being built on a UNIX environment.

Building on Windows follows the same logic but some of the commands may be slightly different.

Build instructions
==================

hictk can be compiled on most UNIX-like systems (including many Linux distributions, MacOS) and Windows.

Build requirements
------------------

Compiling hictk requires a compiler toolchain supporting C++17, such as:

* GCC 8+
* Clang 8+
* Apple-Clang 10.0+
* MSVC 19.12+

Based on our testing, hictk binaries compiled on Linux using Clang are noticeably faster than those compiled with GCC.
For this reason we recommend building hictk using a modern version of Clang whenever possible.

Furthermore, the following tools are required:

* CMake 3.25+
* Conan 2+
* git 2.7+
* make or ninja
* Python3.6+ (including :code:`pip`, required to install Conan)


We recommend installing CMake and Conan in a Python `virtualenv <https://virtualenvwrapper.readthedocs.io/en/stable/>`_, but you are of course free to install build dependencies in any way you want.

.. code-block:: bash

  python3 -m venv /tmp/venv
  /tmp/venv/bin/python3 -m pip install pip setuptools --upgrade
  /tmp/venv/bin/python3 -m pip install 'cmake>=3.25' 'conan>=2' ninja

  # NOTE: It's important to activate the venv after installing CMake
  . /tmp/venv/bin/activate

  whereis cmake  # cmake: /tmp/venv/bin/cmake
  whereis conan  # conan: /tmp/venv/bin/conan
  whereis ninja  # ninja: /tmp/venv/bin/ninja

  cmake --version
  conan --version

  # Detect compiler toolchain. It is usually a good idea to explicitly set CC and CXX
  # Use CC=gcc CXX=g++ if clang is not installed on your machine
   conan profile detect --force

Getting the source code
-----------------------

Download from the `Release <https://github.com/paulsengroup/hictk/releases>`_ page (recommended).

.. code-block:: bash

  mkdir /tmp/hictk
  curl -L 'https://github.com/paulsengroup/hictk/archive/refs/tags/v2.0.2.tar.gz' | tar --strip-components=1 -C /tmp/hictk -xzf -


Using git.

.. code-block:: bash

  git clone https://github.com/paulsengroup/hictk.git /tmp/hictk

  cd /tmp/hictk
  git checkout v2.0.2  # Skip this step if you want to build the latest commit from main

Compiling hictk
---------------

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  # Set these variables to the number of CPU cores available on your machine
  # You can check this with e.g.
  # python -c 'import multiprocessing as mp; print(mp.cpu_count())'
  export CONAN_CPU_COUNT=8
  export CMAKE_BUILD_PARALLEL_LEVEL=8
  export CMAKE_POLICY_VERSION_MINIMUM=3.5

  # Install/build dependencies with Conan
  conan install --build=missing \
                -pr default \
                -s build_type=Release \
                -s compiler.cppstd=17 \
                --output-folder=./build/ \
                .

  # Do not pass -G Ninja if you want CMake to use make instead of ninja
  # Use clang whenever possible, as that usually leads to significantly faster hictk binaries.
  # If clang is not installed on your machine, then replace clang and clang++ with e.g. gcc and g++
  # -DCMAKE_C_COMPILER=gcc
  # -DCMAKE_CXX_COMPILER=g++
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_PREFIX_PATH="$PWD/build" \
        -DHICTK_ENABLE_TESTING=ON \
        -DHICTK_ENABLE_FUZZY_TESTING=ON \
        -DHICTK_BUILD_TOOLS=ON \
        -DHICTK_BUILD_BENCHMARKS=OFF \
        -DHICTK_BUILD_EXAMPLES=OFF \
        -G Ninja \
        -S /tmp/hictk \
        -B /tmp/hictk/build

  cmake --build /tmp/hictk/build

  # If you are compiling hictk on Windows you need to pass the build config as well
  # cmake --build /tmp/hictk/build --config Release

To override the default compiler used by CMake, pass the following arguments to the first CMake command: :code:`-DCMAKE_C_COMPILER=path/to/cc -DCMAKE_CXX_COMPILER=path/to/c++`

We highly recommend using the same compiler when running Conan and CMake.

Running automated tests
=======================

The steps outlined in this section are optional but highly recommended.

Unit tests
----------

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  cd /tmp/hictk
  ctest --test-dir build/   \
        --schedule-random   \
        --output-on-failure \
        --no-tests=error    \
        --timeout 120       \
        -j8  # Change this to the number of available CPU cores

A successful run of the test suite will produce an output like the following:

.. code-block:: console

  user@dev:/tmp/hictk$ ctest --test-dir build/ ...
  ...
   96/106 Test  #62: Cooler: dataset linear iteration - LONG ...........................   Passed    2.26 sec
   97/106 Test #104: Transformers (hic) - SHORT ........................................   Passed    2.81 sec
   98/106 Test   #7: Balancing: SCALE (gw) - SHORT .....................................   Passed    2.49 sec
   99/106 Test  #17: Balancing: SCALE (edge cases) - MEDIUM ............................   Passed    2.78 sec
  100/106 Test  #15: Balancing: ICE (inter) - MEDIUM ...................................   Passed    3.17 sec
  101/106 Test   #6: Balancing: SCALE (inter) - SHORT ..................................   Passed    3.06 sec
  102/106 Test   #8: Balancing: AtomicBitSet - SHORT ...................................   Passed    3.52 sec
  103/106 Test  #66: Cooler: utils merge - LONG ........................................   Passed    3.88 sec
  104/106 Test  #61: Cooler: dataset random iteration - MEDIUM .........................   Passed   10.41 sec
  105/106 Test  #63: Cooler: dataset large read/write - LONG ...........................   Passed   12.10 sec
  106/106 Test  #92: HiC: HiCFileWriter - LONG .........................................   Passed   13.03 sec
  100% tests passed, 0 tests failed out of 106

  Total Test time (real) = 101.97 sec

**All tests are expected to pass. Do not ignore test failures!**

.. raw:: html

   <details>
   <summary><a>Troubleshooting test failures</a></summary>

If one or more tests fail, try the following troubleshooting steps before reaching out for help.

#. Make sure you are running :code:`ctest` from the root of the source tree (:code:`/tmp/hictk` if you are following the instructions).
#. Make sure you are passing the correct build folder to :code:`--test-dir`. Pass the absolute path if necessary (i.e. :code:`--test-dir=/tmp/hictk/build/` if you are following the instructions).
#. Re-run :code:`ctest` with :code:`-j1`. This can be necessary on machines with very little memory (e.g. less than 2GB).
#. Before running :code:`ctest`, create a temporary folder where your user has read-write permissions and where there are at least 100-200MB of space available.
   Then set variable :code:`TMPDIR` to that folder and re-run `ctest`.
#. Checksum the test dataset located under :code:`test/data/` by running :code:`sha256sum -c checksums.sha256`.
   If the checksumming fails or the folder doesn't exist, download and extract the :code:`.tar.zst` file listed in file :code:`cmake/FetchTestDataset.cmake`. Make sure you run :code:`tar -xf` from the root of the repository (:code:`/tmp/hictk` if you are following the instructions).

Example:

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  cd /tmp/hictk

  # Make sure this is the URL listed in file cmake/FetchTestDataset.cmake
  curl -L 'https://zenodo.org/records/13849053/files/hictk_test_data.tar.zst?download=1' | zstdcat | tar -xf -

  # This should print "OK" if the check is successful
  (cd test/data && sha256sum --quiet -c checksums.sha256 && 2>&1 echo OK)

  mkdir ~/hictk-test-dir  # Remember to delete this folder

  TMPDIR="$HOME/hictk-test-dir"      \
  ctest --test-dir=/tmp/hictk/build/ \
        --schedule-random            \
        --output-on-failure          \
        --no-tests=error             \
        --timeout 600                \
        -j1

  # rm -r ~/hictk-test-dir

If after trying the above steps the tests are still failing, feel free to start `discussion <https://github.com/paulsengroup/hictk/discussions>`_ asking for help.

.. raw:: html

   </details>


Integration tests
-----------------

The integration test suite is implemented in Python, requires 3.11 or newer, and can be installed using pip:


.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  pip install test/integration

  hictk_integration_suite --help

Once installed, the full integration suite can be run as follows:

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  cd /tmp/hictk

  hictk_integration_suite \
     build/src/hictk/hictk \
     test/integration/config.toml \
     --data-dir test/data \
     --threads 8 \
     --result-file results.json

  # To run specific parts of the integration suite, pass e.g. --suites=metadata,validate


Installation
============

Once all tests have passed, :code:`hictk` can be installed as follows:

.. code-block:: console

  # Activate venv
  user@dev:/tmp$ . /tmp/venv/bin/activate

  # Install system-wide (requires root/admin rights)
  user@dev:/tmp$ cmake --install /tmp/hictk/build
  -- Install configuration: "Release"
  -- Installing: /usr/local/bin/hictk
  -- Set non-toolchain portion of runtime path of "/usr/local/bin/hictk" to ""
  -- Installing: /usr/local/share/licenses/hictk/LICENSE
  -- Installing: /usr/local/include/hictk
  ...

  # Alternatively, install to custom path
  user@dev:/tmp$ cmake --install /tmp/hictk/build --prefix "$HOME/.local/"
  -- Install configuration: "Release"
  -- Installing: /home/user/.local/bin/hictk
  -- Set non-toolchain portion of runtime path of "/home/user/.local/bin/hictk" to ""
  -- Installing: /home/user/.local/share/licenses/hictk/LICENSE
  -- Installing: /home/user/.local/include/hictk
  ...

  # Install the hictk binary only (i.e. without the header files required for development)
  user@dev:/tmp$ cmake --install /tmp/hictk/build --component Runtime
  -- Install configuration: "Release"
  -- Installing: /usr/local/bin/hictk
  -- Set non-toolchain portion of runtime path of "/usr/local/bin/hictk" to ""
  -- Installing: /usr/local/share/licenses/hictk/LICENSE


Cleaning build artifacts
========================

After successfully compiling hictk the following folders safely be removed:

* Python virtualenv: :code:`/tmp/venv`
* hictk source tree: :code:`/tmp/hictk`

If you are not using Conan in any other project feel free to also delete Conan's folder :code:`~/.conan2/`.
