<!--
Copyright (C) 2024 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# README.md

Test to ensure hictk can be included in a project through `add_subdirectory()`.

```bash
conan install conanfile.py --output-folder build/ --build=missing
ln -s ../../../ hictk_root
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build -DCMAKE_PREFIX_PATH="$PWD/build"
cmake --build build/

build/hictk_test_add_subdirectory
```
