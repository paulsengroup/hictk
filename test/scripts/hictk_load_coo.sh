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

export function readlink_py

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
resolution=10000
batch_size=999999

export PATH="$PATH:$script_dir"

if [ $status -ne 0 ]; then
  exit $status
fi

if ! check_test_files_exist.sh "$ref_cooler"; then
  exit 1
fi

outdir="$(mktemp -d -t hictk-tmp-XXXXXXXXXX)"
trap 'rm -rf -- "$outdir"' EXIT

"$hictk_bin" dump -t chroms "$ref_cooler" > "$outdir/chrom.sizes"

if [[ "$sorted" == true ]]; then
  "$hictk_bin" dump -t pixels "$ref_cooler" |
    "$hictk_bin" load \
      -f coo \
      --assume-sorted \
      --chunk-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
else
  "$hictk_bin" dump -t pixels "$ref_cooler" |
     shuffle.sh |
    "$hictk_bin" load \
      -f coo \
      --assume-unsorted \
      --chunk-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.cool"
fi

if ! compare_matrix_files.sh "$hictk_bin" "$outdir/out.cool" "$ref_cooler" "$resolution"; then
  status=1
fi


if [[ "$sorted" == false ]]; then
  "$hictk_bin" dump -t pixels "$ref_cooler" |
    shuffle.sh |
    "$hictk_bin" load \
      -f coo \
      --assume-unsorted \
      --chunk-size "$batch_size" \
      --bin-size "$resolution" \
      --tmpdir "$outdir" \
      "$outdir/chrom.sizes" \
      "$outdir/out.hic"

  if ! compare_matrix_files.sh "$hictk_bin" "$outdir/out.hic" "$ref_cooler" "$resolution"; then
    status=1
  fi
fi

if [ "$status" -eq 0 ]; then
  printf '\n### PASS ###\n'
else
  printf '\n### FAIL ###\n'
fi

exit "$status"
