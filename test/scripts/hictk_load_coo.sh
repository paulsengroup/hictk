#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

echo "##########################"
echo "#### hictk load (COO) ####"

# readlink -f is not available on macos...
function readlink_py {
  set -eu
  python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$1"
}

function check_files_exist {
  set -eu
  status=0
  for f in "$@"; do
    if [ ! -f "$f" ]; then
      2>&1 echo "Unable to find test file \"$f\""
      status=1
    fi
  done

  return "$status"
}

function compare_coolers {
  set -o pipefail
  set -e

  2>&1 echo "Comparing $1 with $2..."
  if diff <(cooler dump -t chroms "$1") \
          <(cooler dump -t chroms "$2") \
     && \
     diff <(cooler dump --join "$1") \
          <(cooler dump --join "$2");
  then
    2>&1 echo "Files are identical"
    return 0
  else
    2>&1 echo "Files differ"
    return 1
  fi
}

function shuffle {
  if command -v shuf &> /dev/null; then
    shuf
  else
    sort -R
  fi
}

export function readlink_py shuffle

status=0

if [ $# -ne 2 ]; then
  2>&1 echo "Usage: $0 path_to_hictk [un]sorted"
  status=1
fi

hictk_bin="$1"
if [[ "$2" == 'sorted' ]]; then
  sorted=true
else
  sorted=false
fi

data_dir="$(readlink_py "$(dirname "$0")/../data/integration_tests")"
script_dir="$(readlink_py "$(dirname "$0")")"

ref_cooler="$data_dir/4DNFIKNWM36K.subset.fixed-bins.cool"

export PATH="$PATH:$script_dir"

if ! command -v cooler &> /dev/null; then
  2>&1 echo "Unable to find cooler in your PATH"
  status=1
fi

# Try to detect the error outlined below as early as possible:
# https://github.com/open2c/cooler/pull/298
cooler --help > /dev/null

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_files_exist "$ref_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

cooler dump -t chroms "$ref_cooler" > "$outdir/chrom.sizes"

if [[ "$sorted" == true ]]; then
  cooler dump -t pixels "$ref_cooler" |
    "$hictk_bin" load \
      -f coo \
      --assume-sorted \
      --batch-size 1000000 \
      --bin-size 10000 \
      --zero-based \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
else
  cooler dump -t pixels "$ref_cooler" |
     shuffle |
    "$hictk_bin" load \
      -f coo \
      --assume-unsorted \
      --batch-size 1000000 \
      --bin-size 10000 \
      --zero-based \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
fi

if ! compare_coolers "$outdir/out.cool" "$ref_cooler"; then
  status=1
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
