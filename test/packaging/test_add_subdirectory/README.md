# README.md

Test to ensure hictk can be included in a project through `add_subdirectory()`.

```bash
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build -DCMAKE_PREFIX_PATH=$(readlink -f hictk_root/conan-envs/gcc/Release)
cmake --build build/

build/hictk_test_add_subdirectory
```
