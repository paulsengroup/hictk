<!--
Copyright (C) 2024 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictk_fuzzer

The instructions assume `hictk_fuzzer` is being built and run on a Linux machine.
The fuzzer should be capable of running on other OS as well, but some of the steps may require some tweaking.

Dependencies are installed using micromamba, but using PIP and/or your system package manager will do just fine.

## Installing build and runtime dependencies

```bash
micromamba create -c conda-forge -c bioconda -f scripts/environment.yml

# If you don't have a modern C++ compiler toolchain (e.g. GCC8+, Clang 8+, or Apple-clang 11+), you will need to install one.
# Ideally, you should do this with your system package manager.
# If that's not an option, you can try
# micromamba install -c conda-forge -n hictk_fuzzer c-compiler cxx-compiler
```

## Building the fuzzer

```bash
micromamba activate hictk_fuzzer

conan install conanfile.py                  \
   --build='missing'                        \
   -s build_type=Release                    \
   -s compiler.libcxx=libstdc++11           \
   -s compiler.cppstd=17                    \
   --output-folder conan-env

cmake -DCMAKE_BUILD_TYPE=Release           \
      -DCMAKE_PREFIX_PATH="$PWD/conan-env" \
      -DENABLE_DEVELOPER_MODE=ON           \
      -DOPT_ENABLE_CLANG_TIDY=OFF          \
      -DOPT_ENABLE_CPPCHECK=OFF            \
      -DHICTK_ENABLE_TESTING=ON            \
      -DHICTK_ENABLE_FUZZY_TESTING=ON      \
      -DHICTK_BUILD_EXAMPLES=OFF           \
      -DHICTK_DOWNLOAD_TEST_DATASET=OFF    \
      -S ../../                            \
      -B build

cmake --build build -j $(nproc) -t hictk_fuzzer

build/test/fuzzer/src/hictk_fuzzer --help
```

## Running the fuzzer

Running the fuzzer requires a pair of matched files (see later steps).
The reference file should always be a .\[m]cool file.
The test matrix can be in any format supported by hictk.

Furthermore, the fuzzer depends on Cooler

### Installing runtime dependencies

### Fuzzing hictk

```bash
micromamba activate hictk_fuzzer

build/test/fuzzer/src/hictk_fuzzer fuzz \
  --resolution 50000 \
  --nproc $(nproc) \
  --format df \
  test.hic \
  reference.cool
```

## Generating matched datasets from scratch

**IMPORTANT!!** This is usually a slow and demanding process (both in terms of storage, time, and memory).

### Downloading the requirements

Generating matched datasets from scratch requires 4 input files:

- A file with the reference genome in `.chrom.sizes` format (e.g. [hg38.chrom.sizes](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes))
- A file with interactions in `.pairs` format (e.g. [4DNFIYECESRC.pairs.gz](https://data.4dnucleome.org/files-processed/4DNFIYECESRC/#details))
- HiCTools' .jar file: [v3.30.00](https://github.com/aidenlab/HiCTools/releases/tag/v3.30.00)
- JuicerTools' .jar file: [v2.20.00](https://github.com/aidenlab/Juicebox/releases/tag/v2.20.00)

### Generating matched datasets with fixed bin width

```bash
micromamba activate hictk_fuzzer

# Sort/filter .chrom.sizes file

grep '^chr[0-9XY]\+[[:space:]]' hg38.chrom.sizes |
  sort -k1,1V > hg38.filtered.chrom.sizes

# Depending on the size of the pairs file and the requested resolutions/normalizations this step could take several days
scripts/create_datasets_for_fuzzer.py \
  4DNFIYECESRC.pairs.gz \
  4DNFIYECESRC \
  --chrom-sizes hg38.filtered.chrom.sizes \
  --juicer-tools-jar juicer_tools.2.20.00.jar \
  --hic-tools-jar hic_tools.3.30.00.jar \
  --tmpdir ./tmp/ \
  --nproc $(nproc) \
  -Xmx 500G \
  --resolutions     1000 \
                    2000 \
                    5000 \
                   10000 \
                   25000 \
                   50000 \
                  100000 \
                  250000 \
                  500000 \
                 1000000 \
                 2500000 \
                 5000000 \
                10000000
```

### Generating marched datasets with variable bin width

```bash
micromamba activate hictk_fuzzer

# Sort/filter .chrom.sizes file

grep '^chr[0-9XY]\+[[:space:]]' hg38.chrom.sizes |
  sort -k1,1V > hg38.filtered.chrom.sizes

scripts/create_variable_bins_cooler_for_fuzzer.py \
  4DNFIYECESRC.pairs.gz \
  4DNFIYECESRC \
  --bin-size-avg 100000 \
  --bin-size-std 50000 \
  --chrom-sizes hg38.filtered.chrom.sizes \
  --tmpdir ./tmp/ \
  --nproc $(nproc)
```

## Downloading matched datasets used by hictk's CI

Datasets used by hictk's CI are hosted on Zenodo at the following DOIs:

- [10.5281/zenodo.13754754](https://doi.org/10.5281/zenodo.13754754) - .cool datasets
- [10.5281/zenodo.13754785](https://doi.org/10.5281/zenodo.13754785) - .hic v8
- [10.5281/zenodo.13754807](https://doi.org/10.5281/zenodo.13754807) - .hic v9 datasets

Downloading datasets can be automated using `scripts/download_test_datasets.py`

```bash
scripts/download_test_datasets.py \
  test_files.json \
  . \
  --format cool hic8 hic9 \
  --resolution 50000 \
  --dataset 4DNFIYECESRC \
  --nproc 3
```
