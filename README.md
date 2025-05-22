<!--
Copyright (C) 2023 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# hictk

---

<!-- markdownlint-disable MD033 -->

<table>
    <tr>
      <td>Downloads</td>
      <td>
        <a href="https://anaconda.org/bioconda/hictk">
          <img src="https://img.shields.io/conda/vn/bioconda/hictk?label=bioconda&logo=Anaconda" alt="Bioconda">
        </a>
        <a href="https://conan.io/center/recipes/hictk">
          <img src="https://img.shields.io/conan/v/hictk" alt="Conan Center Index">
        </a>
        <a href="https://hub.docker.com/r/paulsengroup/hictk">
          <img src="https://img.shields.io/docker/pulls/paulsengroup/hictk" alt="DockerHub">
        </a>
        <a href="https://doi.org/10.5281/zenodo.8214220">
          <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.8214220.svg" alt="Zenodo">
        </a>
      </td>
    </tr>
    <tr>
      <td>Documentation</td>
      <td>
        <a href="https://hictk.readthedocs.io/">
          <img src="https://readthedocs.org/projects/hictk/badge/?version=latest" alt="Documentation">
        </a>
      </td>
    </tr>
    <tr>
      <td>License</td>
      <td>
        <a href="https://github.com/paulsengroup/hictk/blob/main/LICENSE">
          <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
        </a>
      </td>
    </tr>
    <tr>
      <td>Coverage</td>
      <td>
        <a href="https://codecov.io/gh/paulsengroup/hictk">
          <img src="https://codecov.io/gh/paulsengroup/hictk/graph/badge.svg" alt="Coverage">
        </a>
      </td>
    </tr>
    <tr>
      <td>CI</td>
      <td>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/ubuntu-ci.yml/badge.svg" alt="Ubuntu CI Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/macos-ci.yml/badge.svg" alt="macOS CI Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/windows-ci.yml/badge.svg" alt="Windows CI Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/build-dockerfile.yml/badge.svg" alt="Build Dockerfile Status">
        </a>
      </td>
    </tr>
    <tr>
      <td>CodeQL</td>
      <td>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/codeql-cpp.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/codeql-cpp.yml/badge.svg" alt="CodeQL (C++) Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/codeql-python.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/codeql-python.yml/badge.svg" alt="CodeQL (Python) Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/codeql-actions.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/codeql-actions.yml/badge.svg" alt="CodeQL (GH Actions) Status">
        </a>
      </td>
    </tr>
    <tr>
      <td>Fuzzy Testing</td>
      <td>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/fuzzy-testing.yml/badge.svg" alt="Fuzzy Testing Status">
        </a>
      </td>
    </tr>
    <tr>
      <td>Static Analysis</td>
      <td>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/run-clang-tidy.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/run-clang-tidy.yml/badge.svg" alt="clang-tidy Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/lint-cmakelists.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/lint-cmakelists.yml/badge.svg" alt="Lint CMakeLists.txt files Status">
        </a>
        <a href="https://github.com/paulsengroup/hictk/actions/workflows/lint-cff.yml">
          <img src="https://github.com/paulsengroup/hictk/actions/workflows/lint-cff.yml/badge.svg" alt="Lint CITATION.cff Status">
        </a>
      </td>
    </tr>
</table>

<!-- markdownlint-enable MD033 -->

---

hictk is a blazing fast toolkit to work with .hic and .cool files.

The toolkit consists of a native CLI application and a C++ library running on Linux, macOS, and Windows.\
hictk offers native IO support for Cooler and .hic files, meaning that its implementation is independent of that of cooler, JuicerTools, or straw.

hictk can also be accessed from several programming languages using one of the following libraries:

- [hictkpy](https://github.com/paulsengroup/hictkpy) - Python bindings for hictk: read and write .cool and .hic files directly from Python
- [hictkR](https://github.com/paulsengroup/hictkR) - R bindings for hictk: read .cool and .hic files directly from R
- [libhictk](https://github.com/paulsengroup/hictk) - The native C++ library that underlies hictk

## Features

### Supported formats

The CLI application and C++ library are capable of reading and writing files in the following formats:

<!-- markdownlint-disable MD033 -->

| Format | Revision   | Read | Write           |
| ------ | ---------- | ---- | --------------- |
| .cool  | v1-3 (all) | ✅   | ✅ <sup>1</sup> |
| .mcool | v1-2 (all) | ✅   | ✅ <sup>2</sup> |
| .scool | v1 (all)   | ✅   | ✅ <sup>3</sup> |
| .hic   | v6-9       | ✅   | ✅ <sup>4</sup> |

<small><small>

<sup>1</sup> v3 only\
<sup>2</sup> v2 only\
<sup>3</sup> libhictk only\
<sup>4</sup> v9 only

</small></small>

<!-- markdownlint-enable MD033 -->

### Supported operations

- Seamless conversion between Cooler and .hic formats (from hic to cool and vice versa)
- Uniform interface to query interaction matrices
- High performance and low memory requirements (see benchmarks in the [Supplementary Text](https://academic.oup.com/bioinformatics/article/40/7/btae408/7698028#468964761) from our paper)
- Easy access to file metadata
- Create files from interaction pairs or pre-binned interaction counts
  (e.g. [4DN-DCIC pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md#example-pairs-file)
  or [BEDPE/bedGraph2](https://cooler.readthedocs.io/en/latest/datamodel.html#genomically-labeled-arrays))
- Merge interactions from multiple files into a single file (also supports merging files in different formats)
- Detect (and when possible fix) corrupted files
- Balance interaction matrices using ICE, SCALE, or VC
- Create multi-resolution files suitable for visualization with JuiceBox and HiGlass

All the above operations can be performed on both Cooler and .hic files and yield identical results.

## Installation

hictk can be installed using containers, bioconda, Conan, or directly from source.\
Refer to the [Installation](https://hictk.readthedocs.io/en/stable/installation.html) section in the documentation for more information.

## Quickstart

### hictk (CLI)

hictk provides the following subcommands:

| Subcommand             | Description                                                                                    |
| ---------------------- | ---------------------------------------------------------------------------------------------- |
| **balance**            | Balance Hi-C files using ICE, SCALE, or VC.                                                    |
| **convert**            | Convert Hi-C files between different formats.                                                  |
| **dump**               | Read interactions and other kinds of data from .hic and Cooler files and write them to stdout. |
| **fix-mcool**          | Fix corrupted .mcool files.                                                                    |
| **load**               | Build .cool and .hic files from interactions in various text formats.                          |
| **merge**              | Merge multiple Cooler or .hic files into a single file.                                        |
| **metadata**           | Print file metadata to stdout.                                                                 |
| **rename-chromosomes** | Rename chromosomes found in a Cooler file.                                                     |
| **validate**           | Validate .hic and Cooler files.                                                                |
| **zoomify**            | Convert single-resolution Cooler and .hic files to multi-resolution by coarsening.             |

Refer to the
[Quickstart (CLI)](https://hictk.readthedocs.io/en/stable/quickstart_cli.html) and
[CLI Reference](https://hictk.readthedocs.io/en/stable/cli_reference.html)
sections in the documentation for more details.

### libhictk

libhictk can be installed in various ways, including with Conan and CMake FetchContent.\
Section [Quickstart (API)](https://hictk.readthedocs.io/en/stable/quickstart_api.html) of hictk documentation contains further details on how this can be accomplished.

[Quickstart (API)](https://hictk.readthedocs.io/en/stable/quickstart_api.html) also demonstrates the basic functionality offered by libhictk.\
For more complex examples refer to the sample programs under the [examples/](./examples/) folder as well as to the [source code](./src/hictk/) of hictk.

The public C++ API of hictk is documented in the [C++ API Reference](https://hictk.readthedocs.io/en/stable/cpp_api/index.html) section of hictk documentation.

## Citing

If you use hictk or any of its language bindings in your research, please cite the following publication:

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
