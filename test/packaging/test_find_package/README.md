# README.md

Test to ensure hictk can be included in a project through `find_package()`.

```bash
# Build and install hictk
conan install conanfile.py --output-folder build/
cmake -DCMAKE_BUILD_TYPE=Release -S ../../../ -B hictk_build -DHICTK_ENABLE_TESTING=OFF -DHICTK_BUILD_TOOLS=OFF -DCMAKE_PREFIX_PATH="$PWD/build" -DCMAKE_INSTALL_PREFIX=hictk_install
cmake --build hictk_build/
cmake --install hictk_build/

# Build test project
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build -DCMAKE_PREFIX_PATH="$PWD/build;$PWD/hictk_install/lib/cmake/hictk/"
cmake --build build/

build/hictk_test_find_package
```
