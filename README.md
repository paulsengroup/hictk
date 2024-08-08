<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictk

[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
[![docs](https://readthedocs.org/projects/hictk/badge/?version=latest)](https://hictk.readthedocs.io/en/latest/?badge=latest)
[![Ubuntu CI](https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml)
[![MacOS CI](https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml)
[![Windows CI](https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml)
[![Build Dockerfile](https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml)
[![Fuzzy testing](https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml/badge.svg)](https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml)
[![Download from Bioconda](https://img.shields.io/conda/vn/bioconda/hictk?label=bioconda&logo=Anaconda)](https://anaconda.org/bioconda/hictk)

[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8214220.svg)](https://doi.org/10.5281/zenodo.8214220)

---

hictk is a blazing fast toolkit to work with .hic and .cool files.

This repository hosts `hictk`: a set of CLI tools to work with Cooler, as well as `libhictk`: the C++ library underlying `hictk`.

Python bindings for `libhictk` are available at [paulsengroup/hictkpy](https://github.com/paulsengroup/hictkpy).

hictk is capable of reading files in `.cool`, `.mcool`, `.scool` and `.hic` format (including hic v9) as well as writing `.hic`, `.cool` and `.mcool` files.

## Installing hictk

hictk is developed on Linux and tested on Linux, MacOS and Windows.

hictk can be installed using containers, bioconda or directly from source. Refer to [Installation](https://hictk.readthedocs.io/en/latest/installation.html) for more information.

## Running hictk

hictk provides the following subcommands:

| subcommand             | description                                                                        |
|------------------------|------------------------------------------------------------------------------------|
| __balance__            | Balance HiC matrices using ICE, SCALE or VC.                                       |
| __convert__            | Convert matrices to a different format.                                            |
| __dump__               | Dump data from .hic and Cooler files to stdout.                                    |
| __fix-mcool__          | Fix corrupted .mcool files.                                                        |
| __load__               | Build .cool and .hic files from interactions in various text formats.              |
| __merge__              | Merge multiple Cooler or .hic files into a single file.                            |
| __rename-chromosomes__ | Rename chromosomes found in a Cooler file.                                         |
| __validate__           | Validate .hic and Cooler files.                                                    |
| __zoomify__            | Convert single-resolution Cooler and .hic files to multi-resolution by coarsening. |

Refer to [Quickstart (CLI)](https://hictk.readthedocs.io/en/latest/quickstart_cli.html) and [CLI Reference](https://hictk.readthedocs.io/en/latest/cli_reference.html) for more details.

## Using libhictk

libhictk can be installed in various way, including with Conan and CMake FetchContent. Section [Quickstart (API)](https://hictk.readthedocs.io/en/latest/quickstart_api.html) of hictk documentation contains further details on how this can be accomplished.

[Quickstart (API)](https://hictk.readthedocs.io/en/latest/quickstart_api.html) also showcases the basic functionality offered by libhictk. For more complex examples refer to the sample programs under the [examples/](./examples/) folder as well as to the [source code](./src/hictk/) of hictk.

The public C++ API of hictk is documented in the [C++ API Reference](https://hictk.readthedocs.io/en/latest/cpp_api/index.html) section of hictk documentation.

## Citing

If you use hictk in you reaserch, please cite the following publication:

Roberto Rossini, Jonas Paulsen, hictk: blazing fast toolkit to work with .hic and .cool files
_Bioinformatics_, Volume 40, Issue 7, July 2024, btae408, [https://doi.org/10.1093/bioinformatics/btae408](https://doi.org/10.1093/bioinformatics/btae408)

<details>
<summary>BibTex</summary>

```bibtex
@article{hictk,
    author = {Rossini, Roberto and Paulsen, Jonas},
    title = "{hictk: blazing fast toolkit to work with .hic and .cool files}",
    journal = {Bioinformatics},
    volume = {40},
    number = {7},
    pages = {btae408},
    year = {2024},
    month = {06},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btae408},
    url = {https://doi.org/10.1093/bioinformatics/btae408},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/7/btae408/58385157/btae408.pdf},
}
```

</details>
