#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini (roberros@uio.no)
# SPDX-License-Identifier: MIT

set -e
set -u

# shellcheck disable=SC2064
trap "cd '$PWD'" EXIT

git_root="$(readlink -f "$(git rev-parse --show-toplevel)")"

wd="$git_root/conan-envs"
conanfile="$git_root/conanfile.py"

for compiler in gcc clang; do
  for build_type in Debug Release RelWithDebInfo; do
    CC="$compiler"
    if [[ "$compiler" == gcc* ]]; then
      CXX="${compiler/gcc/g++}${compiler#gcc}"
      profile=gcc
    else
      CXX="${compiler/clang/clang++}${compiler#clang}"
      profile=clang
    fi

    export CC
    export CXX

    outdir="$wd/$compiler/$build_type"
    rm -rf "$outdir"
    mkdir -p "$outdir"

    conan install \
      --build=missing \
      --update \
      -pr "$profile"  \
      -s compiler.cppstd=17 \
      -s build_type="$build_type" \
      --output-folder="$outdir" \
      "$conanfile"

     conan install \
       --build=missing \
       --update \
       -pr "$profile"  \
       -s compiler.cppstd=17 \
       -s build_type="$build_type" \
       -o shared=True \
       --output-folder="$outdir" \
       "$conanfile"
  done
done
