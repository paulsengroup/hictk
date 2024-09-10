<!--
Copyright (C) 2024 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# README.md

Test to ensure hictk can be included in a project through `add_subdirectory()`.

```bash
ln -s ../../../ hictk_root
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build -DCMAKE_PREFIX_PATH=$(readlink -f hictk_root/conan-envs/gcc/Release)
cmake --build build/

build/hictk_test_add_subdirectory
```
