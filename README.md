<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictk

[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
[![Ubuntu CI](https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml)
[![MacOS CI](https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml)
[![Windows CI](https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml)
[![Build Dockerfile](https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml)
[![Fuzzy testing](https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml)
[![Download from Bioconda](https://img.shields.io/conda/vn/bioconda/hictk?label=bioconda&logo=Anaconda)](https://anaconda.org/bioconda/hictk)

<!-- [![Zenodo DOI]()]() -->
---

hictk is a blazing fast toolkit to work with .hic and .cool files.

This repository hosts `hictk`: a set of CLI tools to work with Cooler, as well as `libhictk`: the C++ library underlying `hictk`.

Python bindings for `libhictk` are available at [paulsengroup/hictkpy](https://github.com/paulsengroup/hictkpy).

hictk is capable of reading files in `.cool`, `.mcool`, `.scool` and `.hic` format (including hic v9) as well as writing `.cool` and `.mcool` files.

## Installing hictk

hictk is developed on Linux and tested on Linux, MacOS and Windows.

hictk can be installed or compiled with one of the following methods.

### Docker or Singularity/Apptainer

First, make sure you follow the instructions on how to install Docker or Singularity/Apptainer on your OS.

<details>
<summary>Installing Docker</summary>

The following instructions assume you have root/admin permissions.

- [Linux](https://docs.docker.com/desktop/install/linux-install/#generic-installation-steps/)
- [MacOS](https://docs.docker.com/desktop/install/mac-install/)
- [Windows](https://docs.docker.com/desktop/install/windows-install/)

On some Linux distributions just installing Docker is not enough.
You also need to start (and optionally enable) the appropriate service(s).
This is usually done with one of the following:

```bash
sudo systemctl start docker
sudo systemctl start docker.service
```

Refer to [Docker](https://docs.docker.com/engine/install/) or your distribution documentation for more details.

</details>

<details>
<summary>Installing Singularity/Apptainer</summary>

The following instructions assume you have root/admin permissions.

Apptainer can be easily installed using your system package manager.

[Here](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) you can find instructions for common Linux distributions such as Ubuntu.

Even if your distribution is not listed in the above documentation, your system package manager likely includes a package for Singularity or Apptainer. If this is not the case, then you must install Apptainer from source (instructions available [here](https://github.com/apptainer/apptainer/blob/release-1.1/INSTALL.md)).

</details>

#### Pulling hictk Docker image

hictk Docker images are available on [ghcr.io](https://github.com/paulsengroup/hictk/pkgs/container/hictk)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/hictk).

Downloading and running the latest stable release can be done as follows:

```console
# Using Docker, may require sudo
user@dev:/tmp$ docker run ghcr.io/paulsengroup/hictk:0.0.2 --help

# Using Singularity/Apptainer
user@dev:/tmp$ singularity run ghcr.io/paulsengroup/hictk:0.0.2 --help

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
```

The above will print hictk's help message, and is equivalent to running `hictk --help` on the command line (assuming hictk is available on your machine).


### Conda (bioconda)

hictk package for Linux, MacOS and Windows is available on [bioconda](https://anaconda.org/bioconda/hictk) and can be installed as follows:

```console
user@dev:/tmp$ conda create -n hictk -c conda-forge -c bioconda hictk

(hictk) user@dev:/tmp$ conda activate hictk

(hictk) user@dev:/tmp$ whereis hictk
hictk: /home/user/.miniconda3/envs/hictk/bin/hictk

(hictk) user@dev:/tmp$ hictk --version
hictk-v0.0.2-bioconda
```

### Installing from source


hictk can be compiled on most UNIX-like systems, including many Linux distributions, MacOS (10.15+) and Windows.

<details>
<summary>Build instructions</summary>

## Build instructions

Instructions assume hictk is being built on a UNIX environment.

Building on Windows follows the same logic but some of the commands may be slightly different.

### Build requirements

Compiling hictk requires a compiler toolchain supporting C++17, such as:

- GCC 8+
- Clang 8+
- Apple-Clang 10.0+

Furthermore, the following tools are required:
- CMake 3.25+
- Conan 2+
- git 2.7+
- make or ninja
- Python3.6+ (including `pip`, required to install Conan)


We recommend installing CMake and Conan in a Python [virtualenv](https://virtualenvwrapper.readthedocs.io/en/stable/), but you are of course free to install the build dependencies in any way you want.

```bash
python3 -m venv /tmp/venv
/tmp/venv/bin/python3 -m pip install pip setuptools --upgrade
/tmp/venv/bin/python3 -m pip  install 'cmake>=3.25' 'conan>=2' ninja

# NOTE: It's important to activate the venv after installing CMake
. /tmp/venv/bin/activate

whereis cmake  # cmake: /tmp/venv/bin/cmake
whereis conan  # conan: /tmp/venv/bin/conan
whereis ninja  # ninja: /tmp/venv/bin/ninja

cmake --version
conan --version

# Detect compiler toolchain. It is usually a good idea to explicitly set CC and CXX
CC=gcc CXX=g++ conan profile detect --force
```

#### Getting the source code

Download from the [Release](https://github.com/paulsengroup/hictk/releases) page (recommended).
```bash
mkdir /tmp/hictk
curl -L 'https://github.com/paulsengroup/hictk/archive/refs/tags/v0.0.2.tar.gz' | tar --strip-components=1 -C /tmp/hictk -xzf -
```

Using git.
```bash
git clone https://github.com/paulsengroup/hictk.git /tmp/hictk

cd /tmp/hictk
git checkout v0.0.2  # Skip this step if you want to build the latest commit from main
```

#### Compiling hictk

```bash
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
```

To override the default compiler used by CMake, pass the following arguments to the first CMake command: `-DCMAKE_C_COMPILER=path/to/cc -DCMAKE_CXX_COMPILER=path/to/c++`

We highly recommend using the same compiler when running Conan and CMake.

## Running automated tests

The steps outlined in this section are optional but highly recommended.

#### Unit tests

```bash
# Activate venv
. /tmp/venv/bin/activate

cd /tmp/hictk
ctest --test-dir build/   \
      --schedule-random   \
      --output-on-failure \
      --no-tests=error    \
      --timeout 120       \
      -j8  # Change this to the number of available CPU cores
```

A successful run of the test suite will produce an output like the following:
```console
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
```

__All tests are expected to pass. Do not ignore test failures!__

<details>
<summary> Troubleshooting test failures </summary>
If one or more tests fail, try the following troubleshooting steps before reaching out for help.

1. Make sure you are running `ctest` from the root of the source tree (`/tmp/hictk` if you are following the instructions).
2. Make sure you are passing the correct build folder to `--test-dir`. Pass the absolute path if necessary (i.e. `--test-dir=/tmp/hictk/build/` if you are following the instructions).
3. Re-run `ctest` with `-j1`. This can be necessary on machines with very little memory (e.g. less than 2GB).
4. Before running `ctest`, create a temporary folder where your user has read-write permissions and where there are at least 100-200MB of space available.
   Then set variable `TMPDIR` to that folder and re-run `ctest`.
5.  Checksum the test dataset located under `test/data/` by running `sha256sum -c checksums.sha256`.
    If the checksumming fails or the folder doesn't exist, download and extract the `.tar.xz` file listed in file `cmake/FetchTestDataset.cmake`. Make sure you run `tar -xf` from the root of the repository (`/tmp/hictk` if you are following the instructions).

Example:
```bash
# Activate venv
. /tmp/venv/bin/activate

cd /tmp/hictk

# Make sure this is the URL listed  in file cmake/FetchTestDataset.cmake
curl -L 'https://zenodo.org/record/8143316/files/hictk_test_data.tar.xz?download=1' | tar -xJf -

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
```

If after trying the above steps the tests are still failing, feel free to start [discussion](https://github.com/paulsengroup/hictk/discussions) asking for help.

</details>


#### Integration tests

The integration test scripts depend on the following tools:

- cooler>=0.9
- java
- [juicer_tools](https://github.com/aidenlab/Juicebox/releases/latest) or [hic_tools](https://github.com/aidenlab/HiCTools/releases/latest)
- xz
- common UNIX shell commands

cooler can be installed using pip:
```bash
/tmp/venv/bin/pip3 install 'cooler>=0.9'
```

juicer_tools and hic_tools do not need to be installed, downloading the JAR file is enough:
```bash
curl -L 'https://github.com/aidenlab/HiCTools/releases/download/v3.30.00/hic_tools.3.30.00.jar' -o /tmp/hictk/hic_tools.jar
```

If not already installed, `xz` can usually be installed with your system package manager (on some Linux distributions the relevant package is called `xz-utils`).

```bash
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
```

## Installation

Once all tests have passed, `hictk` can be installed as follows:

```console
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
```

## Cleaning build artifacts

After successfully compiling hictk the following folders safely be removed:
- Python virtualenv: `/tmp/venv`
- hictk source tree: `/tmp/hictk`

If you are not using Conan in any other project feel free to also delete Conan's folder `~/.conan2/`

</details>


### Running hictk

hictk provides the following subcommands:

| subcommand   | description                                                                            |
|--------------|----------------------------------------------------------------------------------------|
| __convert__  | Convert matrices between .hic and Cooler formats.                                      |
| __dump__     | Write interactions from .hic or Cooler files to the terminal.                          |
| __load__     | Generate a Cooler file from pixels or pairs of interactions in text format.            |
| __merge__    | Merge multiple Cooler files using the same reference assembly.                         |
| __validate__ | Validate Cooler and .hic files.                                                        |
| __zoomify__  | Convert single-resolution cooler files to multi-resolution cooler files by coarsening. |

#### Examples

Converting hic->cooler:

```bash
# Create a .mcool file using all resolutions available in interactions.hic
hictk convert interactions.hic interactions.mcool

# Create a .cool file at 10kb resolution
hictk convert interactions.hic interactions.cool --resolutions 10000

# Create a .mcool file using a subset of the resolutions available in interactions.hic
hictk convert interactions.hic interactions.mcool --resolutions 10000 20000 50000
```

Converting cool->hic:

```bash
# Create a .hic file using interactions.cool as base resolution
hictk convert interactions.cool interactions.hic --juicer-tools-jar /tmp/hic_tools.jar

# Create a .hic file with the resolutions found in interactions.mcool
hictk convert interactions.mcool interactions.cool --juicer-tools-jar /tmp/hic_tools.jar
```

Dumping interactions:

```shell
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
```

Loading interactions in a Cooler file:

```bash
# Create a 10kbp .cool file using hg38 as reference
hictk load --format 4dn --assembly hg38 hg38.chrom.sizes 10000 out.cool < interactions.txt

# Same as above but using gzip-compressed interactions
zcat interactions.txt.gz | hictk load --format 4dn --assembly hg38 hg38.chrom.sizes 10000 out.cool

# Using interactions in bedgraph2 format (see --help for the list of supported formats)
hictk load --format bg2 --assembly hg38 hg38.chrom.sizes 10000 out.cool < interactions.txt
```

Merging multiple coolers:

```bash
hictk merge interactions1.cool interactions2.cool -o merged.cool
```

Checking file integrity (especially useful to detect corrupted .mcool from 4DNucleome, see [here](https://github.com/robomics/20221129_4dnucleome_bug_report) and [here](https://github.com/open2c/cooler/issues/319)):

```bash
hictk validate interactions.cool --validate-index

hictk validate interactions.hic
```

Creating .mcool files from .cool files:

```bash
hictk zoomify interactions.cool interactions.mcool --resolutions 1000 5000 10000 ...

# Coarsen a single resolution
hictk zoomify interactions.cool interactions.ccool --no-copy-base-resolution --resolutions 10000
```
