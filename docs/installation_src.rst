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
  CC=gcc CXX=g++ conan profile detect --force

Getting the source code
-----------------------

Download from the `Release <https://github.com/paulsengroup/hictk/releases>`_ page (recommended).

.. code-block:: bash

  mkdir /tmp/hictk
  curl -L 'https://github.com/paulsengroup/hictk/archive/refs/tags/v0.0.6.tar.gz' | tar --strip-components=1 -C /tmp/hictk -xzf -


Using git.

.. code-block:: bash

  git clone https://github.com/paulsengroup/hictk.git /tmp/hictk

  cd /tmp/hictk
  git checkout v0.0.12  # Skip this step if you want to build the latest commit from main

Compiling hictk
---------------

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  # Set these variables to the number of CPU cores available on your machine
  # You can check this with e.g.
  # python -c 'import multiprocessing as mp; print(mp.cpu_count())')
  export CONAN_CPU_COUNT=8
  export CMAKE_BUILD_PARALLEL_LEVEL=8

  # Install/build dependencies with Conan
  conan install --build=missing \
                -pr default \
                -s build_type=Release \
                -s compiler.cppstd=17 \
                --output-folder=./build/ \
                .

  # This may take a while, as CMake will run Conan to build hictk dependencies.
  # Do not pass -G Ninja if you want CMake to use make instead of ninja
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH="$PWD/build" \
        -DHICTK_ENABLE_TESTING=ON \
        -DHICTK_BUILD_TOOLS=ON \
        -G Ninja \
        -S /tmp/hictk \
        -B /tmp/hictk/build

  cmake --build /tmp/hictk/build

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
  63/70 Test #21: Cooler: init files - SHORT .......................................   Passed    0.02 sec
  64/70 Test #57: HiC: pixel selector fetch (observed NONE BP 10000) - LONG ........   Passed    1.53 sec
  65/70 Test  #5: Cooler: index validation - SHORT .................................   Passed    3.83 sec
  66/70 Test #17: Cooler: index validation - SHORT .................................   Passed    3.62 sec
  67/70 Test #37: Cooler: utils merge - LONG .......................................   Passed    4.35 sec
  68/70 Test #67: Transformers (cooler) - SHORT ....................................   Passed    4.11 sec
  69/70 Test #36: Cooler: dataset random iteration - MEDIUM ........................   Passed    5.50 sec
  70/70 Test #40: Cooler: dataset large read/write - LONG ..........................   Passed   11.47 sec

  100% tests passed, 0 tests failed out of 70

  Total Test time (real) =  12.03 sec

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
   If the checksumming fails or the folder doesn't exist, download and extract the :code:`.tar.xz` file listed in file :code:`cmake/FetchTestDataset.cmake`. Make sure you run :code:`tar -xf` from the root of the repository (:code:`/tmp/hictk` if you are following the instructions).

Example:

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  cd /tmp/hictk

  # Make sure this is the URL listed in file cmake/FetchTestDataset.cmake
  curl -L 'https://zenodo.org/records/10522583/files/hictk_test_data.tar.xz?download=1' | tar -xJf -

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

The integration test scripts depend on the following tools:

* cooler>=0.9
* java
* `juicer_tools <https://github.com/aidenlab/Juicebox/releases/latest>`_ or `hic_tools <https://github.com/aidenlab/HiCTools/releases/latest>`_
* xz
* common UNIX shell commands

cooler can be installed using pip:

.. code-block:: bash

  /tmp/venv/bin/pip3 install 'cooler>=0.9'

juicer_tools and hic_tools do not need to be installed, downloading the JAR file is enough:

.. code-block:: bash

  curl -L 'https://github.com/aidenlab/HiCTools/releases/download/v3.30.00/hic_tools.3.30.00.jar' -o /tmp/hictk/hic_tools.jar

If not already installed, :code:`xz` can usually be installed with your system package manager (on some Linux distributions the relevant package is called :code:`xz-utils`).

.. code-block:: bash

  # Activate venv
  . /tmp/venv/bin/activate

  cd /tmp/hictk

  # hictk convert
  test/scripts/hictk_convert_cool2hic.sh build/src/hictk/hictk juicer_tools.jar
  test/scripts/hictk_convert_hic2cool.sh build/src/hictk/hictk

  # hictk dump
  test/scripts/hictk_dump_balanced.sh build/src/hictk/hictk
  test/scripts/hictk_dump_bins.sh build/src/hictk/hictk
  test/scripts/hictk_dump_chroms.sh build/src/hictk/hictk
  test/scripts/hictk_dump_cis.sh build/src/hictk/hictk
  test/scripts/hictk_dump_gw.sh build/src/hictk/hictk
  test/scripts/hictk_dump_trans.sh build/src/hictk/hictk

  # hictk load (sorted)
  test/scripts/hictk_load_4dn.sh build/src/hictk/hictk sorted
  test/scripts/hictk_load_bg2.sh build/src/hictk/hictk sorted
  test/scripts/hictk_load_coo.sh build/src/hictk/hictk sorted

  # hictk load (unsorted)
  test/scripts/hictk_load_4dn.sh build/src/hictk/hictk unsorted
  test/scripts/hictk_load_bg2.sh build/src/hictk/hictk unsorted
  test/scripts/hictk_load_coo.sh build/src/hictk/hictk unsorted

  # hictk merge
  test/scripts/hictk_merge.sh build/src/hictk/hictk

  # hictk validate
  test/scripts/hictk_validate.sh build/src/hictk/hictk

  # hictk zoomify
  test/scripts/hictk_zoomify.sh build/src/hictk/hictk

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
  -- Set runtime path of "/usr/local/bin/hictk" to ""
  -- Up-to-date: /usr/local/share/licenses/hictk/LICENSE
  ...

  # Alternatively, install to custom path
  user@dev:/tmp$ cmake --install /tmp/hictk/build --prefix "$HOME/.local/"
  -- Install configuration: "Release"
  -- Installing: /home/user/.local/bin/hictk
  -- Set runtime path of "/home/user/.local/bin/hictk" to ""
  -- Up-to-date: /home/user/.local/share/licenses/hictk/LICENSE
  ...

Cleaning build artifacts
========================

After successfully compiling hictk the following folders safely be removed:

* Python virtualenv: :code:`/tmp/venv`
* hictk source tree: :code:`/tmp/hictk`

If you are not using Conan in any other project feel free to also delete Conan's folder :code:`~/.conan2/`
